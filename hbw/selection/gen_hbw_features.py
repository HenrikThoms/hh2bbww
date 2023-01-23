# coding: utf-8

"""
Selectors to set ak columns for gen particles of hh2bbww
"""

from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, Route, EMPTY_FLOAT
from columnflow.selection import Selector, selector
from columnflow.selection import SelectionResult
from columnflow.production import Producer, producer

ak = maybe_import("awkward")

coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


def masked_sorted_indices(mask: ak.Array, sort_var: ak.Array, ascending: bool = False) -> ak.Array:
    """
    Helper function to obtain the correct indices of an object mask
    """
    indices = ak.argsort(sort_var, axis=-1, ascending=ascending)
    return indices[mask[indices]]


@selector(
    uses={
        "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.btagDeepFlavB",
        "gen_hbw_decay",
    },
    produces=set(
        f"cutflow.{gp}_{var}"
        for gp in ["h1", "h2", "b1", "b2", "wlep", "whad", "l", "nu", "q1", "q2", "sec1", "sec2",]
        for var in ["pt", "eta", "phi", "mass"]
    ) | set(
        f"cutflow.{var}_{p1}_{p2}"
        for var in ["dR", "dETA", "dPHI"]
        for p1, p2 in [("sec1", "sec2"), ("b1", "b2"), ("h1", "h2"), ("wlep", "whad"), ("q1", "q2"), ("l", "nu")]
    ) | {
        "cutflow.dETArel_sec1_sec2", "cutflow.m_sec1_sec2", "cutflow.m_h1_h2",
    } | {
        "cf_matched_m_fjfj", "cf_matched_dR_fjfj", "cf_matched_dEta_fjfj",
    } | {
        "cf_matched_forward_jet1_pt", "cf_matched_forward_jet2_pt",
    }
)
def gen_hbw_decay_features(self: Selector, events: ak.Array, **kwargs) -> ak.Array:

    # enable 4-vector behavior
    for gp in ["h1", "h2", "b1", "b2", "wlep", "whad", "l", "nu", "q1", "q2", "sec1", "sec2"]:
        events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
        events["gen_hbw_decay", gp] = ak.with_name(events.gen_hbw_decay[gp], "PtEtaPhiMLorentzVector")

    for var in ["pt", "eta", "phi", "mass"]:
        for gp in ["h1", "h2", "b1", "b2", "wlep", "whad", "l", "nu", "q1", "q2", "sec1", "sec2",]:
            events = set_ak_column(events, f"cutflow.{gp}_{var}", events.gen_hbw_decay[gp][var])
    # Forwardjet matching
    forward_jet_mask = (events.Jet.pt > 5) & (abs(events.Jet.eta < 4.7))
    fjet_indices = masked_sorted_indices(forward_jet_mask, events.Jet.pt)
    fjets = events.Jet[fjet_indices]

    sec1 = events.gen_hbw_decay["sec1"]
    sec2 = events.gen_hbw_decay["sec2"]
    
    dr1 = fjets.delta_r(sec1)
    dr2 = fjets.delta_r(sec2)

    fjets["dr_sec1"] = dr1
    fjets["dr_sec2"] = dr2

    sec1_match = (dr1 < 0.4)
    sec2_match = (dr2 < 0.4)

    multiple_matched_sec1 = fjets[sec1_match]
    multiple_matched_sec2 = fjets[sec2_match]

    matched_sec1 = ak.pad_none(multiple_matched_sec1[ak.argsort(multiple_matched_sec1.dr_sec1)], 1)[:, 0]
    matched_sec2 = ak.pad_none(multiple_matched_sec2[ak.argsort(multiple_matched_sec2.dr_sec2)], 1)[:, 0]

    #events = set_ak_column(events, "cf_matched_sec1", matched_sec1)
    """
    config.add_variable(
        name="cf_matched_sec1_pt",
        expression="cf_matched_sec1.pt",
    )
    """
    events = set_ak_column(events, "matched_sec1", matched_sec1)
    events = set_ak_column(events, "matched_sec2", matched_sec2)

    events = set_ak_column(events, f"cf_matched_forward_jet1_pt", matched_sec1.pt)
    events = set_ak_column(events, f"cf_matched_forward_jet2_pt", matched_sec2.pt)

    m_matched_fjfj =(matched_sec1 + matched_sec2).mass
    events = set_ak_column(events, "cf_matched_m_fjfj", ak.fill_none(m_matched_fjfj, EMPTY_FLOAT))

    dR_matched_fjfj = matched_sec1.delta_r(matched_sec2)
    events = set_ak_column(events, "cf_matched_dR_fjfj", ak.fill_none(dR_matched_fjfj, EMPTY_FLOAT))

    dEta_matched_fjfj = abs(matched_sec1.eta - (matched_sec2.eta))
    events = set_ak_column(events, "cf_matched_dEta_fjfj", ak.fill_none(dEta_matched_fjfj, EMPTY_FLOAT))

    # 4-vector example

    for p1, p2 in [("sec1", "sec2"), ("b1", "b2"), ("h1", "h2"), ("wlep", "whad"), ("q1", "q2"), ("l", "nu")]:
        var = events.gen_hbw_decay[p1].delta_r(events.gen_hbw_decay[p2])
        events = set_ak_column(events, f"cutflow.dR_{p1}_{p2}", var)
        var2 = events.gen_hbw_decay[p1].delta_phi(events.gen_hbw_decay[p2])
        events = set_ak_column(events, f"cutflow.dPHI_{p1}_{p2}", var2)
        var3 = abs(events.gen_hbw_decay[p1].eta - events.gen_hbw_decay[p2].eta)
        events = set_ak_column(events, f"cutflow.dETA_{p1}_{p2}", var3)

    dETArel_sec1_sec2 = (abs(events.gen_hbw_decay.sec1.eta) -abs(events.gen_hbw_decay.sec2.eta))

    events = set_ak_column(events, "cutflow.dETArel_sec1_sec2", dETArel_sec1_sec2)
        
    m_sec1_sec2 = (events.gen_hbw_decay.sec1 + events.gen_hbw_decay.sec2).mass

    events = set_ak_column(events, "cutflow.m_sec1_sec2", m_sec1_sec2)

    m_h1_h2 = (events.gen_hbw_decay.h1 + events.gen_hbw_decay.h2).mass

    events = set_ak_column(events, "cutflow.m_h1_h2", m_h1_h2)

    #dR_sec1_sec2 = events.gen_hbw_decay.sec1.delta_r(events.gen_hbw_decay.sec2)

    #events = set_ak_column(events, "cutflow.dR_sec1_sec2", dR_sec1_sec2)

    #dETA_sec1_sec2 = abs(events.gen_hbw_decay.sec1.eta - events.gen_hbw_decay.sec2.eta)

    #events = set_ak_column(events, "cutflow.dETA_sec1_sec2", dETA_sec1_sec2)

    #dPHI_sec1_sec2 = events.gen_hbw_decay.sec1.delta_phi(events.gen_hbw_decay.sec2)

    #events = set_ak_column(events, "cutflow.dPHI_sec1_sec2", dPHI_sec1_sec2)

    return events


@selector(uses={"event",  "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.btagDeepFlavB",
                "gen_hbw_decay", gen_hbw_decay_features})
def sel_fjet1(self: Selector, events: ak.Array, results: SelectionResult, **kwargs) -> ak.Array:
    # self[gen_hbw_decay_features](events)
    return(~ak.is_none(events.matched_sec1))

@selector(uses={"event",  "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.btagDeepFlavB",
                "gen_hbw_decay", gen_hbw_decay_features})
def sel_fjet2(self: Selector, events: ak.Array, results: SelectionResult, **kwargs) ->ak.Array:
    return(~ak.is_none(events.matched_sec2))


@producer(
    uses={
        "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.genJetIdx",
        "GenJet.pt", "GenJet.eta", "GenJet.phi", "GenJet.mass", "GenJet.partonFlavour", "GenJet.hadronFlavour",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass",
        "gen_hbw_decay",
    },
)
def gen_hbw_matching(
        self: Producer, events: ak.Array,
        results: SelectionResult = None, verbose: bool = False, 
        dR_req: float = 0.4, ptdiff_req: float = 10.,
        **kwargs,
) -> ak.Array:
    """
    Function that matches HH->bbWW decay product gen particles to Reco-level jets and leptons.
    """

    gen_matches = {}

    # jet matching
    # NOTE: might be nice to remove jets that already have been matched
    for gp_tag in ("b1", "b2", "q1", "q2", "sec1", "sec2"):
        gp = events.gen_hbw_decay[gp_tag]

        dr = events.Jet.delta_r(gp)
        pt_diff = (events.Jet.pt - gp.pt) / gp.pt

        jet_match_mask = (dr < dR_req) & (abs(pt_diff) < ptdiff_req)
        jet_matches = events.Jet[jet_match_mask]

        if verbose:
            print(gp_tag, "multiple matches:", ak.sum(ak.num(jet_matches) > 1))
            print(gp_tag, "no matches:", ak.sum(ak.num(jet_matches) == 0))

        # if multiple matches found: choose match with smallest dr
        jet_match = jet_matches[ak.argsort(jet_matches.delta_r(gp))]
        gen_matches[gp_tag] = jet_match

    # lepton matching for combined electron and muon
    lepton_fields = set(events.Muon.fields).intersection(events.Electron.fields)

    lepton = ak.concatenate([
        ak.zip({f: events.Muon[f] for f in lepton_fields}),
        ak.zip({f: events.Electron[f] for f in lepton_fields}),
    ], axis=-1)
    lepton = ak.with_name(lepton, "PtEtaPhiMLorentzVector")

    gp_lep = events.gen_hbw_decay.l
    dr = lepton.delta_r(gp_lep)
    pt_diff = (lepton.pt - gp_lep.pt) / gp_lep.pt
    lep_match_mask = (dr < dR_req) & (abs(pt_diff) < ptdiff_req)
    lep_matches = lepton[lep_match_mask]

    if verbose:
        print("l multiple matches:", ak.sum(ak.num(lep_matches) > 1))
        print("l no matches:", ak.sum(ak.num(lep_matches) == 0))

    lep_match = lep_matches[ak.argsort(lep_matches.delta_r(gp_lep))]
    gen_matches["l"] = lep_match

    # write matches into events and return them
    events = set_ak_column(events, "gen_match", ak.zip(gen_matches, depth_limit=1))

    return events
