# coding: utf-8

"""
Definition of variables.
"""

import order as od

from columnflow.columnar_util import EMPTY_FLOAT


def add_variables(config: od.Config) -> None:
    """
    Adds all variables to a *config*.
    """
    config.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(200, -10, 10),
        x_title="MC weight",
    )

    # Event properties
    config.add_variable(
        name="n_jet",
        binning=(11, -0.5, 10.5),
        x_title="Number of jets",
    )
    config.add_variable(
        name="n_deepjet",
        binning=(11, -0.5, 10.5),
        x_title="Number of deepjets",
    )
    config.add_variable(
        name="n_electron",
        binning=(11, -0.5, 10.5),
        x_title="Number of electrons",
    )
    config.add_variable(
        name="n_muon",
        binning=(11, -0.5, 10.5),
        x_title="Number of muons",
    )
    config.add_variable(
        name="ht",
        binning=(40, 0, 1500),
        x_title="HT",
    )
    config.add_variable(
        name="ht_rebin",
        expression="ht",
        binning=[0, 80, 120, 160, 200, 240, 280, 320, 400, 500, 600, 800, 1200],
        unit="GeV",
        x_title="HT",
    )

    # Object properties

    # Jets (4 pt-leading jets)
    for i in range(4):
        config.add_variable(
            name=f"jet{i+1}_pt",
            expression=f"Jet.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0., 400.),
            unit="GeV",
            x_title=r"Jet %i $p_{T}$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_eta",
            expression=f"Jet.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(50, 0., 5),
            x_title=r"Jet %i $\eta$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_phi",
            expression=f"Jet.phi[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, -3.2, 3.2),
            x_title=r"Jet %i $\phi$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_mass",
            expression=f"Jet.mass[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, -3.2, 3.2),
            x_title=r"Jet %i mass" % (i + 1),
        )

    # Bjets (2 b-score leading jets)
#    for i in range(2):
#        config.add_variable(
#            name=f"bjet{i+1}_pt",
#            expression=f"Bjet.pt[:,{i}]",
#            null_value=EMPTY_FLOAT,
#            binning=(40, 0., 400.),
#            unit="GeV",
#            x_title=r"Bjet %i $p_{T}$" % (i + 1),
#        )
#        config.add_variable(
#            name=f"bjet{i+1}_eta",
#            expression=f"Bjet.eta[:,{i}]",
#            null_value=EMPTY_FLOAT,
#            binning=(50, 0., 5),
#            x_title=r"Bjet %i $\eta$" % (i + 1),
#        )
#        config.add_variable(
#            name=f"bjet{i+1}_phi",
#            expression=f"Jet.phi[:,{i}]",
#            null_value=EMPTY_FLOAT,
#            binning=(40, -3.2, 3.2),
#            x_title=r"Bjet %i $\phi$" % (i + 1),
#        )
#        config.add_variable(
#            name=f"bjet{i+1}_mass",
#            expression=f"Bjet.mass[:,{i}]",
#            null_value=EMPTY_FLOAT,
#            binning=(40, -3.2, 3.2),
#            x_title=r"Bjet %i mass" % (i + 1),
#        )

    # Leptons
    for obj in ["Electron", "Muon"]:
        config.add_variable(
            name=f"{obj.lower()}_pt",
            expression=f"{obj}.pt[:,0]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0., 400.),
            unit="GeV",
            x_title=obj + r" $p_{T}$",
        )
        config.add_variable(
            name=f"{obj.lower()}_phi",
            expression=f"{obj}.phi[:,0]",
            null_value=EMPTY_FLOAT,
            binning=(40, -3.2, 3.2),
            x_title=obj + r" $\phi$",
        )
        config.add_variable(
            name=f"{obj.lower()}_eta",
            expression=f"{obj}.eta[:,0]",
            null_value=EMPTY_FLOAT,
            binning=(50, 0., 5),
            x_title=obj + r" $\eta$",
        )
        config.add_variable(
            name=f"{obj.lower()}_mass",
            expression=f"{obj}.mass[:,0]",
            null_value=EMPTY_FLOAT,
            binning=(40, -3.2, 3.2),
            x_title=obj + " mass",
        )

    # MET
    config.add_variable(
        name="met_pt",
        expression="MET.pt[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(40, 0., 400.),
        unit="GeV",
        x_title=r"MET $p_{T}$",
    )
    config.add_variable(
        name="met_phi",
        expression="MET.phi[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(40, -3.2, 3.2),
        x_title=r"MET $\phi$",
    )

    # bb features
    config.add_variable(
        name="m_bb",
        binning=(40, 0., 400.),
        unit="GeV",
        x_title=r"$m_{bb}$",
    )
    config.add_variable(
        name="deltaR_bb",
        binning=(40, 0, 5),
        x_title=r"$\Delta R(b,b)$",
    )
    # jj features
    config.add_variable(
        name="m_jj",
        binning=(40, 0., 400.),
        unit="GeV",
        x_title=r"$m_{jj}$",
    )
    config.add_variable(
        name="deltaR_jj",
        binning=(40, 0, 5),
        x_title=r"$\Delta R(j_{1},j_{2})$",
    )
##################################################
    # Gen particles

    # cutflow variables
    for i in range(6):
        config.add_variable(
            name=f"cf_jet{i+1}_BtagFlavB",
            expression=f"cf_jet{i+1}_btag",
            binning=(40, 0, 1),
            unit="GeV",
            x_title=r"Btag Flavour Score of Jets %s" % (i+1),
        )
        config.add_variable(
            name=f"cf_jet{i+1}_pt",
            expression=f"cf_jet{i+1}_pt",
            binning=(40, 0, 400),
            unit="GeV",
            x_title=r"$p_{T}$ of Jet %s" % (i+1),
        )
    config.add_variable(
        name="cf_invariante_m_fjfj",
        expression="cf_m_fjfj",
        binning=(40, 0, 300),
        x_title=r"Invariant mass of the forward Jets",
    )
    config.add_variable(
        name="cf_2dinvariante_m_fjfj",
        expression="cf_m_fjfj",
        binning=(40, 0, 240),
        x_title=r"Invariant mass of the forward Jets",
    )
    config.add_variable(
        name="cf_dR_fjfj",
        expression="cf_dR_fjfj",
        binning=(40, 0, 12),
        x_title=r"$\Delta$R between the forward Jets",
    )
    config.add_variable(
        name="cf_dEta_fjfj",
        expression="cf_dEta_fjfj",
        binning=(40, 0, 9),
        x_title=r"$\Delta\eta$ between the forward Jets",
    )
    config.add_variable(
        name="cf_2ddEta_fjfj",
        expression="cf_dEta_fjfj",
        binning=(40, 0, 3),
        x_title=r"$\Delta\eta$ between the forward Jets",
    )
    config.add_variable(
        name="cf_dPhi_fjfj",
        expression="cf_dPhi_fjfj",
        binning=(40, -3.5, 3.5),
        x_title=r"$\Delta\phi$ between the forward Jets",
    )
    for i in range(2):
        config.add_variable(
            name=f"cf_forward_jet{i+1}_pt",
            expression=f"cf_forward_jet{i+1}_pt",
            binning=(40, 0, 250),
            x_title=r"Pt of Forward Jet %s" % (i+1), 
        )
        config.add_variable(
            name=f"cf_forward_jet{i+1}_eta",
            expression=f"cf_forward_jet{i+1}_eta",
            binning=(40, -6, 6),
            x_title=r"$\eta$ of Forward Jet %s" % (i+1), 
        )
        config.add_variable(
            name=f"cf_forward_jet{i+1}_phi",
            expression=f"cf_forward_jet{i+1}_phi",
            binning=(40, -3.5, 3.5),
            x_title=r"$\phi$ of Forward Jet %s" % (i+1), 
        )
        config.add_variable(
            name=f"cf_forward_jet{i+1}_mass",
            expression=f"cf_forward_jet{i+1}_mass",
            binning=(40, 0, 10),
            x_title=r"mass of Forward Jet %s" % (i+1), 
        )
    config.add_variable(
        name="cf_n_jet",
        expression="cf_n_jet",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of central jets",
    )
    config.add_variable(
        name="cf_n_Bjet",
        expression="cf_n_Bjet",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of B jets",
    )
    config.add_variable(
        name="cf_n_electron",
        expression="cf_n_electron",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of electrons",
    )
    config.add_variable(
        name="cf_n_muon",
        expression="cf_n_muon",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of muons",
    )
    config.add_variable(
        name="cf_n_forwardjet",
        expression="cf_n_forwardjet",
        binning=(11, -0.5, 10.5),
        x_title=r"Number of forward jets",
    )
    for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        config.add_variable(
            name=f"cf_dEta_{i}_{j}",
            expression=f"cf_dEta_{i}_{j}",
            binning=(40, 0, 8),
            x_title=r"$\Delta$ $\eta$ of jet %s and %s" % (i+1, j+1),
        )
    for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        config.add_variable(
            name=f"cf_2ddEta_{i}_{j}",
            expression=f"cf_dEta_{i}_{j}",
            binning=(40, 0, 4),
            x_title=r"$\Delta$ $\eta$ of jet %s and %s" % (i+1, j+1),
        )
  #  config.add_variable(
  #      name="cf_dEta_all",
  #      expression="cf_dEta_all",
  #      binning=(40, 0, 12),
  #      x_title=r"$\Delta$ $\eta$ of all combined jets"
  #  )
    for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        config.add_variable(
            name=f"cf_m_{i}_{j}",
            expression=f"cf_m_{i}_{j}",
            binning=(40, 0, 400),
            x_title=r"invariant mass of jet %s and %s" % (i+1, j+1),
        )
  #  config.add_variable(
  #      name="cf_m_all",
  #      expression="cf_m_all",
  #      binning=(40, 0, 10),
  #      x_title=r"invariante mass of all combined jets"
  #  )

    for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        config.add_variable(
            name=f"cf_dPhi_{i}_{j}",
            expression=f"cf_dPhi_{i}_{j}",
            binning=(40, -3.5, 3.5),
            x_title=r"$\Delta$ $\phi$ of jet %s and %s" % (i+1, j+1),
        )
#################################################################################
    for gp in ["h1", "h2", "b1", "b2", "wlep", "whad", "l", "nu", "q1", "q2", "sec1", "sec2"]:
        config.add_variable(
            name=f"gen_{gp}_pt",
            expression=f"cutflow.{gp}_pt",
            binning=(40, 0, 250),
            unit="GeV",
            x_title=r"$p_{T, %s}^{gen}$" % (gp),
        )
        mass_bins = (40, 0, 400)
        if gp in ["b1", "b2"]:
            mass_bins = (40, 0, 1)
        elif gp in ["h1","h2"]:
            mass_bins = (40, 122, 128)
        elif gp in ["l"]:
            mass_bins = (40, 0, 1)
        elif gp in ["nu"]:
            mass_bins = (40, 0, 1)
        elif gp in ["q1", "q2"]:
            mass_bins = (40, 0, 1)
        elif gp in ["sec1", "sec2"]:
            mass_bins = (40, 0, 1)
        elif gp in ["wlep", "whad"]:
            mass_bins = (40, 0, 100)
        else:
            mass_bins = (40, 0, 400)

        config.add_variable(
            name=f"gen_{gp}_mass",
            expression=f"cutflow.{gp}_mass",
            binning=mass_bins,
            unit="GeV",
            x_title=r"$m_{%s}^{gen}$" % (gp),
        )
        config.add_variable(
            name=f"gen_{gp}_eta",
            expression=f"cutflow.{gp}_eta",
            binning=(40, -6., 6.),
            x_title=r"$\eta_{%s}^{gen}$" % (gp),
        )
        config.add_variable(
            name=f"gen_{gp}_phi",
            expression=f"cutflow.{gp}_phi",
            binning=(40, -3.5, 3.5),
            x_title=r"$\phi_{%s}^{gen}$" % (gp),
        )

    for p1, p2 in [("sec1", "sec2"), ("b1", "b2"), ("h1", "h2"), ("wlep", "whad"), ("q1", "q2"), ("l", "nu")]:
    # example dR variable
        config.add_variable(
        name=f"gen_dR_{p1}_{p2}",
        expression=f"cutflow.dR_{p1}_{p2}",
        binning=(40, 0, 12),
        x_title=r"$\Delta$R between %s and %s" % (p1, p2),
        )
        config.add_variable(
            name=f"gen_dETA_{p1}_{p2}",
            expression=f"cutflow.dETA_{p1}_{p2}",
            binning=(40, 0, 12),
            x_title=r"$\Delta\eta$ between %s and %s" % (p1, p2),
        )
        config.add_variable(
            name=f"gen_dPHI_{p1}_{p2}",
            expression=f"cutflow.dPHI_{p1}_{p2}",
            binning=(40, 0, 3.5),
            x_title=r"$\Delta\phi$ between %s and %s" % (p1, p2),
        )

    config.add_variable(
        name="gen_dETArel_sec1_sec2",
        expression="cutflow.dETArel_sec1_sec2",
        binning=(40, 0, 12),
        x_title=r"relative $\eta$ between initial-state partons",
    )

    config.add_variable(
        name="gen_m_sec1_sec2",
        expression="cutflow.m_sec1_sec2",
        binning=(40, 0, 10),
        unit="GeV",
        x_title=r"invariant mass between initial-state partons",
    )
    config.add_variable(
        name="gen_m_h1_h2",
        expression="cutflow.m_h1_h2",
        binning=(40, 250, 400),
        unit="GeV",
        x_title=r"invariant mass between higgs-bosons",
    )

##################################

    config.add_variable(
        name="cf_matched_dR_fjfj",
        expression="cf_matched_dR_fjfj",
        binning=(40, 0, 12),
        x_title=r"$\Delta$R between the matched forward Jets",
    )
    config.add_variable(
        name="cf_matched_dEta_fjfj",
        expression="cf_matched_dEta_fjfj",
        binning=(40, 0, 9),
        x_title=r"$\Delta\eta$ between the matched forward Jets",
    )
    config.add_variable(
        name="cf_matched_m_fjfj",
        expression="cf_matched_m_fjfj",
        binning=(40, 0, 500),
        x_title=r"invariant mass between the matched forward Jets",
    )
    config.add_variable(
        name="cf_matched_forward_jet1_pt",
        expression="cf_matched_forward_jet1_pt",
        binning=(40, 0, 200),
        x_title=r"Pt of the matched primary forward Jet",
    )
    config.add_variable(
        name="cf_matched_forward_jet2_pt",
        expression="cf_matched_forward_jet2_pt",
        binning=(40, 0, 200),
        x_title=r"Pt of the matched secondary forward Jet",
    )
