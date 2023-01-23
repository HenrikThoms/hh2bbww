# coding: utf-8

"""
Selection methods for HHtobbWW.
"""

from collections import defaultdict
from typing import Tuple

from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column 
from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.production.categories import category_ids
from columnflow.production.processes import process_ids
from hbw.production.gen_hbw_decay import gen_hbw_decay_products
from hbw.selection.general import jet_energy_shifts, increment_stats
from hbw.selection.cutflow_features import cutflow_features
from hbw.selection.gen_hbw_features import gen_hbw_decay_features, gen_hbw_matching

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(uses={"event"})
def sel_incl(self: Selector, events: ak.Array, **kwargs) -> ak.Array:
    # select all
    return ak.ones_like(events.event)


def masked_sorted_indices(mask: ak.Array, sort_var: ak.Array, ascending: bool = False) -> ak.Array:
    indices = ak.argsort(sort_var, axis=-1, ascending=ascending)
    return indices[mask[indices]]


@selector(
    uses={"Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.btagDeepFlavB"},
    produces={"cutflow.n_jet", "cutflow.n_deepjet_med"},
    shifts={jet_energy_shifts},
    exposed=True,
)
def forward_jet_selection(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> Tuple[ak.Array, SelectionResult]:

    forward_jet_mask = (events.Jet.pt > 30) & (abs(events.Jet.eta < 4.7))
    fjet_indices = masked_sorted_indices(forward_jet_mask, events.Jet.pt)
    fjets = events.Jet[fjet_indices]

    fjet_pairs = ak.combinations(fjets, 2)

    f0 = fjet_pairs[:, :, "0"]
    f1 = fjet_pairs[:, :, "1"]

    fjet_pairs["deta"] = abs(f0.eta - f1.eta)
    fjet_pairs["invmass"] = (f0 + f1).mass

    fjet_mask = (fjet_pairs.deta > 3) & (fjet_pairs.invmass > 500)

    fjet_selection = ak.sum(fjet_mask >= 1, axis=-1) >= 1

    fjet_pairs = masked_sorted_indices

    # build and return selection results plus new columns
    return events, SelectionResult(
        steps={"ForwardJetPair": fjet_selection},
    )


@selector(
    uses={"Jet.pt", "Jet.eta", "Jet.btagDeepFlavB"},
    produces={"cutflow.n_jet", "cutflow.n_deepjet_med"},
    shifts={jet_energy_shifts},
    exposed=True,
)
def jet_selection(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> Tuple[ak.Array, SelectionResult]:
    # HH -> bbWW(qqlnu) jet selection
    # - require at least 3 jets with pt>30, eta<2.4
    # - require at least 1 jet with pt>30, eta<2.4, b-score>0.3040 (Medium WP)

    # jets
    jet_mask_loose = (events.Jet.pt > 5) & abs(events.Jet.eta < 2.4)
    jet_mask = (events.Jet.pt > 30) & (abs(events.Jet.eta) < 2.4)
    events = set_ak_column(events, "cutflow.n_jet", ak.sum(jet_mask, axis=1))
    jet_sel = ak.sum(jet_mask, axis=1) >= 3
    jet_indices = masked_sorted_indices(jet_mask, events.Jet.pt)

    # forward jets
    forward_jet_mask = (events.Jet.pt > 30) & (abs(events.Jet.eta) >= 2.4) & (abs(events.Jet.eta) < 4.7)
    forward_jet_indices = masked_sorted_indices(forward_jet_mask, events.Jet.pt)

    # b-tagged jets, medium working point
    wp_med = self.config_inst.x.btag_working_points.deepcsv.medium
    bjet_mask = (jet_mask) & (events.Jet.btagDeepFlavB >= wp_med)
    events = set_ak_column(events, "cutflow.n_deepjet_med", ak.sum(bjet_mask, axis=1))
    bjet_sel = ak.sum(bjet_mask, axis=1) >= 1

    # sort jets after b-score and define b-jets as the two b-score leading jets
    bjet_indices = masked_sorted_indices(jet_mask, events.Jet.btagDeepFlavB)[:, :2]
    bjetcut_indices = masked_sorted_indices(bjet_mask, events.Jet.btagDeepFlavB)

    # define b-jets as the two b-score leading jets, b-score sorted
    bjet_indices = masked_sorted_indices(jet_mask, events.Jet.btagDeepFlavB)[:, :2]

    # define lightjets as all non b-jets, pt-sorted
    b_idx = ak.fill_none(ak.pad_none(bjet_indices, 2), -1)
    lightjet_indices = jet_indices[(jet_indices != b_idx[:, 0]) & (jet_indices != b_idx[:, 1])]

    # example column: high jet multiplicity region (>=6 jets)
    events = set_ak_column(events, "jet_high_multiplicity", ak.sum(jet_mask, axis=1) >= 6)

    # build and return selection results plus new columns
    return events, SelectionResult(
        steps={"Jet": jet_sel, "Bjet": bjet_sel},
        objects={"Jet": {
            "Jet": jet_indices,
            "ForwardJet": forward_jet_indices, "Bjet": bjet_indices, "Lightjet": lightjet_indices, "Bjetcut": bjetcut_indices, "LooseJet": masked_sorted_indices(jet_mask_loose, events.Jet.pt),
        }},
    )


@selector(uses={
    "Electron.pt", "Electron.eta", "Electron.cutBased",
    "Muon.pt", "Muon.eta", "Muon.tightId", "Muon.looseId", "Muon.pfRelIso04_all",
    "HLT.Ele35_WPTight_Gsf", "HLT.IsoMu27",
})
def lepton_selection(
        self: Selector,
        events: ak.Array,
        stats: defaultdict,
        **kwargs,
) -> Tuple[ak.Array, SelectionResult]:
    # HH -> bbWW(qqlnu) lepton selection
    # - require exactly 1 lepton (e or mu) with pt_e>28 / pt_mu>25, eta<2.4 and tight ID
    # - require that events are triggered by SingleMu or SingleEle trigger

    # 2017 Trigger selection (TODO different triggers based on year of data-taking)
    trigger_sel = (events.HLT.Ele35_WPTight_Gsf) | (events.HLT.IsoMu27)

    # Lepton definition for this analysis
    e_mask = (events.Electron.pt > 36) & (abs(events.Electron.eta) < 2.4) & (events.Electron.cutBased == 4)
    mu_mask = (
        (events.Muon.pt > 28) &
        (abs(events.Muon.eta) < 2.4) &
        (events.Muon.tightId) &
        (events.Muon.pfRelIso04_all < 0.15)
    )

    lep_sel = ak.sum(e_mask, axis=-1) + ak.sum(mu_mask, axis=-1) == 1
    mu_sel = (ak.sum(e_mask, axis=-1) == 0) & (ak.sum(mu_mask, axis=-1) == 1)

    # determine the masked lepton indices
    e_indices = masked_sorted_indices(e_mask, events.Electron.pt)
    mu_indices = masked_sorted_indices(mu_mask, events.Muon.pt)

    # build and return selection results plus new columns
    return events, SelectionResult(
        steps={
            "Lepton": lep_sel, "Trigger": trigger_sel, "Muon": mu_sel,
        },
        objects={"Electron": {"Electron": e_indices}, "Muon": {"Muon": mu_indices}},
    )


@selector(
    uses={
        jet_selection, forward_jet_selection, lepton_selection, cutflow_features,
        category_ids, process_ids, increment_stats, attach_coffea_behavior,
        "mc_weight",  # not opened per default but always required in Cutflow tasks
    },
    produces={
        jet_selection, forward_jet_selection, lepton_selection, cutflow_features,
        category_ids, process_ids, increment_stats, attach_coffea_behavior,
        "mc_weight",  # not opened per default but always required in Cutflow tasks
    },
    shifts={
        jet_energy_shifts,
    },
    exposed=True,
)
def default(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> Tuple[ak.Array, SelectionResult]:
    # ensure coffea behavior
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # jet selection
    events, jet_results = self[jet_selection](events, stats, **kwargs)
    results += jet_results

    # forward-jet selection
    events, forward_jet_results = self[forward_jet_selection](events, stats, **kwargs)
    results += forward_jet_results

    # lepton selection
    events, lepton_results = self[lepton_selection](events, stats, **kwargs)
    results += lepton_results

    # combined event selection after all steps
    event_sel = (
        jet_results.steps.Jet &
        jet_results.steps.Bjet &
        lepton_results.steps.Lepton &
        lepton_results.steps.Trigger
    )
    results.main["event"] = event_sel

    # build categories
    events = self[category_ids](events, results=results, **kwargs)

    # create process ids
    events = self[process_ids](events, **kwargs)

    # add cutflow features
    events = self[cutflow_features](events, results=results, **kwargs)

    # increment stats
    self[increment_stats](events, event_sel, stats, **kwargs)


    return events, results


@selector(
    uses={
        default, attach_coffea_behavior, "mc_weight",  # mc_weight should be included from default
        gen_hbw_decay_products, gen_hbw_decay_features,  gen_hbw_matching,
    },
    produces={
        attach_coffea_behavior,
        category_ids, process_ids, increment_stats, "mc_weight",
        gen_hbw_decay_products, gen_hbw_decay_features, gen_hbw_matching,
    },
    exposed=True,
)
def gen_hbw(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> Tuple[ak.Array, SelectionResult]:
    """
    Selector that is used to perform GenLevel studies but also allow categorization and event selection
    using the default reco-level selection.
    Should only be used for HH samples
    """

    # ensure coffea behavior
    events = self[attach_coffea_behavior](events, **kwargs)

    # extract relevant gen HH decay products
    events = self[gen_hbw_decay_products](events, **kwargs)

    # produce relevant columns
    events = self[gen_hbw_decay_features](events, **kwargs)

    # run the default Selector
    events, results = self[default](events, stats, **kwargs)

    # match genparticles with reco objects
    events = self[gen_hbw_matching](events, results, **kwargs)

    return events, results


@gen_hbw.init
def gen_hbw_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) and not self.dataset_inst.x("is_hbw", False):
        raise Exception("This selector is only usable for HH samples")

    cat_names = [c.name for c in self.config_inst.categories]
    if "fjet1" not in cat_names:
        self.config_inst.add_category(
            name="fjet1",
            id=300,
            selection="sel_fjet1",
            label="first forward Jet matched correctly",
        )
        
    cat_names = [c.name for c in self.config_inst.categories]
    if "fjet2" not in cat_names:
        self.config_inst.add_category(
            name="fjet2",
            id=301,
            selection="sel_fjet2",
            label="second forward Jet matched correctly",
        )
