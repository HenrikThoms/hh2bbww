# coding: utf-8

"""
Selectors to set ak columns for cutflow features
"""

from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, Route, EMPTY_FLOAT
from columnflow.selection import Selector, SelectionResult, selector

ak = maybe_import("awkward")


@selector(
    uses={"Jet.pt", "Jet.phi", "Jet.eta", "Jet.mass", "Jet.btagDeepFlavB",},
    produces={
        "cf_m_fjfj", "cf_dR_fjfj", "cf_dEta_fjfj", "cf_dPhi_fjfj",# "cf_dEta_all", "cf_m_all",
    } | set(
        f"cf_jet{i+1}_{var}"
        for var in ["pt", "btag"]
        for i in range(6)
    ) | set(
        f"cf_forward_jet{i+1}_{var2}"
        for var2 in ["pt", "eta", "phi", "mass"]
        for i in range(2)
    ) | set(
        f"cf_n_{obj}"
        for obj in ["jet", "electron", "muon", "forwardjet", "Bjet", "Bjetcut"]
    ) | set(
        f"cf_dEta_{i}_{j}"
        for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    ) | set(
        f"cf_m_{i}_{j}"
        for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    )
      | set(
        f"cf_dPhi_{i}_{j}"
        for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    )
        
)
def cutflow_features(self: Selector, events: ak.Array, results: SelectionResult, **kwargs) -> ak.Array:

    # determine jet pt before applying jet pt cut (and ideally after applying eta cut?)
    jet_indices = ak.argsort(events.Jet.pt, ascending=False)
    jets = events.Jet[jet_indices]
    for i in range(6):
        events = set_ak_column(events, f"cf_jet{i+1}_pt", Route(f"pt[:, {i}]").apply(jets, EMPTY_FLOAT))
        events = set_ak_column(events, f"cf_jet{i+1}_btag", Route(f"btagDeepFlavB[:, {i}]").apply(jets, EMPTY_FLOAT))

    forward_jet_indices = results.objects.Jet.ForwardJet
    forward_jets = events.Jet[forward_jet_indices]
    forward_jets = ak.pad_none(forward_jets, 2)
    
    for var in ["pt", "eta", "phi", "mass"]:
        for i in range(2):
            events = set_ak_column(events, f"cf_forward_jet{i+1}_{var}", Route(f"{var}[:, {i}]").apply(forward_jets, EMPTY_FLOAT))

    m_fjfj =(forward_jets[:, 0] + forward_jets[:, 1]).mass
    events = set_ak_column(events, "cf_m_fjfj", ak.fill_none(m_fjfj, EMPTY_FLOAT))

    dR_fjfj = forward_jets[:, 0].delta_r(forward_jets[:, 1])
    events = set_ak_column(events, "cf_dR_fjfj", ak.fill_none(dR_fjfj, EMPTY_FLOAT))

    dEta_fjfj = abs(forward_jets[:, 0].eta - (forward_jets[:, 1].eta))
    events = set_ak_column(events, "cf_dEta_fjfj", ak.fill_none(dEta_fjfj, EMPTY_FLOAT))
    
    dPhi_fjfj = forward_jets[:, 0].delta_phi(forward_jets[:, 1])
    events = set_ak_column(events, "cf_dPhi_fjfj", ak.fill_none(dPhi_fjfj, EMPTY_FLOAT))

    # Number of objects should be counted after appyling
    events = set_ak_column(events, "cf_n_jet", ak.num(results.objects.Jet.Jet, axis=1))
    events = set_ak_column(events, "cf_n_electron", ak.num(results.objects.Electron.Electron, axis=1))
    events = set_ak_column(events, "cf_n_muon", ak.num(results.objects.Muon.Muon, axis=1))
    events = set_ak_column(events, "cf_n_forwardjet", ak.num(results.objects.Jet.ForwardJet, axis=1))
    events = set_ak_column(events, "cf_n_Bjet", ak.num(results.objects.Jet.Bjet, axis=1))
    events = set_ak_column(events, "cf_n_Bjetcut", ak.num(results.objects.Jet.Bjetcut, axis=1))
    # delta eta over *all* jet combinations

    jet_pairs = ak.combinations(events.Jet, 2, fields=("first", "second"))
    #dEta_all =abs(jet_pairs.first.eta - jet_pairs.second.eta)
    #m_all = (jet_pairs.first + jet_pairs.second).mass

    # works here, but it's possible the plotting this won't work
    # because there is more than one entry per event

    #events = set_ak_column(events, "cf_dEta_all", ak.fill_none(dEta_all, EMPTY_FLOAT))
    #events = set_ak_column(events, "cf_m_all", ak.fill_none(m_all, EMPTY_FLOAT))

    # helper array to get the index pairs of each combination
    jet_index_pairs = ak.argcombinations(events.Jet, 2, fields=("first", "second"))

    # select certain jet combinations explicitly
    for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        jet_pairs_ij = ak.firsts(
            jet_pairs[
                (jet_index_pairs.first == i) & 
                (jet_index_pairs.second == j)
            ],
            axis=-1
        )

        dEta_ij = abs(jet_pairs_ij.first.eta - jet_pairs_ij.second.eta)
        events = set_ak_column(events, f"cf_dEta_{i}_{j}", ak.fill_none(dEta_ij, EMPTY_FLOAT))

        m_ij = (jet_pairs_ij.first + jet_pairs_ij.second).mass
        events = set_ak_column(events, f"cf_m_{i}_{j}", ak.fill_none(m_ij, EMPTY_FLOAT))

        dPhi_ij = jet_pairs_ij.first.delta_phi(jet_pairs_ij.second)
        events = set_ak_column(events, f"cf_dPhi_{i}_{j}", ak.fill_none(dPhi_ij, EMPTY_FLOAT))
        
   # for var in ["pt", "eta"]:
    #    for i in range(2):
     #       events = set_ak_column(events, f"cf_Muon{i+1}_{var}", Route(f"{var}[:, {i}]").apply(Muon, EMPTY_FLOAT))
            
  #  for var in ["pt", "eta"]:
   #     for i in range(2):
    #        events = set_ak_column(events, f"cf_Electron{i+1}_{var}", Route(f"{var}[:, {i}]").apply(Electron, EMPTY_FLOAT))


    return events
