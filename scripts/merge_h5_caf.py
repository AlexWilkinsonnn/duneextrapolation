"""
Add true and predicted FD reco to the original ND CAF file. Since this will be working with /pnfs
locations, the files are copied to the working directory, processed, and then copied back out
presumably to a /pnfs location.
"""
import subprocess, os, argparse
from array import array

import h5py
import ROOT
import numpy as np

def main(args):
    ret = subprocess.run(["ifdh", "ls", args.h5_dir], stdout=subprocess.PIPE)
    f_h5_paths = ret.stdout.decode("utf-8").strip("\n").split("\n")[1:] # first is directory

    for f_h5_path in f_h5_paths:
        f_h5_name, f_caf_name = gather_inputs(f_h5_path, args.caf_dir)

        try:
            f_h5 = h5py.File(f_h5_name)
        except Exception as e:
            print("Failed to open h5 file with error:")
            print(e)
            print("Skipping")
            gather_outputs(f_caf_name, f_h5_name)
            continue

        if "insidendlar_fd_reco" not in f_h5.keys() or "predresp_fd_reco" not in f_h5.keys():
            print("h5 file missing datasets, skipping")
            gather_outputs(f_caf_name, f_h5_name)
            continue

        f_caf = ROOT.TFile.Open(f_caf_name, "UPDATE")
        t_caf = f_caf.Get("caf")
        t_eventid = f_caf.Get("eventid")
        t_fdreco = ROOT.TTree("FDRecoFriend", "FDRecoFriend")

        # General info
        b_eventid = array("i", [0])
        t_fdreco.Branch("eventID", b_eventid, "eventID/I")

        # Pred CVN scores
        b_pred_cvn_anu = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultAntineutrino", b_pred_cvn_anu, "FDPredCVNResultAntineutrino/F"
        )
        b_pred_cvn_nue = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultNue", b_pred_cvn_nue, "FDPredCVNResultNue/F"
        )
        b_pred_cvn_numu = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultNumu", b_pred_cvn_numu, "FDPredCVNResultNumu/F"
        )
        b_pred_cvn_nc = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultNC", b_pred_cvn_nc, "FDPredCVNResultAntineutrino/F"
        )
        b_pred_cvn_nutau = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultNutau", b_pred_cvn_nutau, "FDPredCVNResultNutau/F"
        )
        b_pred_cvn_0p = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult0p", b_pred_cvn_0p, "FDPredCVNResult0p/F"
        )
        b_pred_cvn_1p = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult1p", b_pred_cvn_1p, "FDPredCVNResult1p/F"
        )
        b_pred_cvn_2p = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult2p", b_pred_cvn_2p, "FDPredCVNResult2p/F"
        )
        b_pred_cvn_np = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultnp", b_pred_cvn_np, "FDPredCVNResultnp/F"
        )
        b_pred_cvn_0pi = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult0pi", b_pred_cvn_0pi, "FDPredCVNResult0pi/F"
        )
        b_pred_cvn_1pi = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult1pi", b_pred_cvn_1pi, "FDPredCVNResult1pi/F"
        )
        b_pred_cvn_2pi = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult2pi", b_pred_cvn_2pi, "FDPredCVNResult2pi/F"
        )
        b_pred_cvn_npi = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultnpi", b_pred_cvn_npi, "FDPredCVNResultnpi/F"
        )
        b_pred_cvn_0pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult0pi0", b_pred_cvn_0pi0, "FDPredCVNResult0pi0/F"
        )
        b_pred_cvn_1pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult1pi0", b_pred_cvn_1pi0, "FDPredCVNResult1pi0/F"
        )
        b_pred_cvn_2pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult2pi0", b_pred_cvn_2pi0, "FDPredCVNResult2pi0/F"
        )
        b_pred_cvn_npi0 = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultnpi0", b_pred_cvn_npi0, "FDPredCVNResultnpi0/F"
        )
        b_pred_cvn_0n = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult0n", b_pred_cvn_0n, "FDPredCVNResult0n/F"
        )
        b_pred_cvn_1n = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult1n", b_pred_cvn_1n, "FDPredCVNResult1n/F"
        )
        b_pred_cvn_2n = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResult2n", b_pred_cvn_2n, "FDPredCVNResult2n/F"
        )
        b_pred_cvn_nn = array("f", [0])
        t_fdreco.Branch(
            "FDPredCVNResultnn", b_pred_cvn_nn, "FDPredCVNResultnn/F"
        )

        # Pred energies
        b_pred_numu_reco_method = array("i", [0])
        t_fdreco.Branch(
            "FDPredEvRecoNumuMethod", b_pred_numu_reco_method, "FDPredEvRecoNumuMethod/I"
        )
        b_pred_numu_nu_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoNumu", b_pred_numu_nu_E, "FDPredEvRecoNumu/F"
        )
        b_pred_numu_lep_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoLepNumu", b_pred_numu_lep_E, "FDPredEvRecoLepNumu/F"
        )
        b_pred_numu_had_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoHadNumu", b_pred_numu_had_E, "FDPredEvRecoHadNumu/F"
        )
        b_pred_nc_nu_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoNC", b_pred_nc_nu_E, "FDPredEvRecoNC/F"
        )
        b_pred_nc_lep_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoLepNC", b_pred_nc_lep_E, "FDPredEvRecoLepNC/F"
        )
        b_pred_nc_had_E = array("f", [0])
        t_fdreco.Branch(
            "FDPredEvRecoHadNC", b_pred_nc_had_E, "FDPredEvRecoHadNC/F"
        )

        # Pred hit info
        b_pred_n_hits_Z = array("i", [0])
        t_fdreco.Branch(
            "FDPredNumHitsZ", b_pred_n_hits_Z, "FDPredNumHitsZ/I"
        )
        b_pred_sum_hits_summedadc_Z = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsSummedADCZ", b_pred_sum_hits_summedadc_Z, "FDPredSumHitsSummedADCZ/F"
        )
        b_pred_sum_hits_integral_Z = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsIntegralZ", b_pred_sum_hits_integral_Z, "FDPredSumHitsIntegralZ/F"
        )
        b_pred_n_hits_U = array("i", [0])
        t_fdreco.Branch(
            "FDPredNumHitsU", b_pred_n_hits_U, "FDPredNumHitsU/I"
        )
        b_pred_sum_hits_summedadc_U = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsSummedADCU", b_pred_sum_hits_summedadc_U, "FDPredSumHitsSummedADCU/F"
        )
        b_pred_sum_hits_integral_U = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsIntegralU", b_pred_sum_hits_integral_U, "FDPredSumHitsIntegralU/F"
        )
        b_pred_n_hits_V = array("i", [0])
        t_fdreco.Branch(
            "FDPredNumHitsV", b_pred_n_hits_V, "FDPredNumHitsV/I"
        )
        b_pred_sum_hits_summedadc_V = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsSummedADCV", b_pred_sum_hits_summedadc_V, "FDPredSumHitsSummedADCV/F"
        )
        b_pred_sum_hits_integral_V = array("f", [0])
        t_fdreco.Branch(
            "FDPredSumHitsIntegralV", b_pred_sum_hits_integral_V, "FDPredSumHitsIntegralV/F"
        )

        # True CVN scores
        b_true_cvn_anu = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultAntineutrino", b_true_cvn_anu, "FDTrueCVNResultAntineutrino/F"
        )
        b_true_cvn_nue = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultNue", b_true_cvn_nue, "FDTrueCVNResultNue/F"
        )
        b_true_cvn_numu = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultNumu", b_true_cvn_numu, "FDTrueCVNResultNumu/F"
        )
        b_true_cvn_nc = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultNC", b_true_cvn_nc, "FDTrueCVNResultAntineutrino/F"
        )
        b_true_cvn_nutau = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultNutau", b_true_cvn_nutau, "FDTrueCVNResultNutau/F"
        )
        b_true_cvn_0p = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult0p", b_true_cvn_0p, "FDTrueCVNResult0p/F"
        )
        b_true_cvn_1p = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult1p", b_true_cvn_1p, "FDTrueCVNResult1p/F"
        )
        b_true_cvn_2p = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult2p", b_true_cvn_2p, "FDTrueCVNResult2p/F"
        )
        b_true_cvn_np = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultnp", b_true_cvn_np, "FDTrueCVNResultnp/F"
        )
        b_true_cvn_0pi = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult0pi", b_true_cvn_0pi, "FDTrueCVNResult0pi/F"
        )
        b_true_cvn_1pi = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult1pi", b_true_cvn_1pi, "FDTrueCVNResult1pi/F"
        )
        b_true_cvn_2pi = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult2pi", b_true_cvn_2pi, "FDTrueCVNResult2pi/F"
        )
        b_true_cvn_npi = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultnpi", b_true_cvn_npi, "FDTrueCVNResultnpi/F"
        )
        b_true_cvn_0pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult0pi0", b_true_cvn_0pi0, "FDTrueCVNResult0pi0/F"
        )
        b_true_cvn_1pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult1pi0", b_true_cvn_1pi0, "FDTrueCVNResult1pi0/F"
        )
        b_true_cvn_2pi0 = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult2pi0", b_true_cvn_2pi0, "FDTrueCVNResult2pi0/F"
        )
        b_true_cvn_npi0 = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultnpi0", b_true_cvn_npi0, "FDTrueCVNResultnpi0/F"
        )
        b_true_cvn_0n = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult0n", b_true_cvn_0n, "FDTrueCVNResult0n/F"
        )
        b_true_cvn_1n = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult1n", b_true_cvn_1n, "FDTrueCVNResult1n/F"
        )
        b_true_cvn_2n = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResult2n", b_true_cvn_2n, "FDTrueCVNResult2n/F"
        )
        b_true_cvn_nn = array("f", [0])
        t_fdreco.Branch(
            "FDTrueCVNResultnn", b_true_cvn_nn, "FDTrueCVNResultnn/F"
        )

        # True energies
        b_true_numu_reco_method = array("i", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoNumuMethod", b_true_numu_reco_method, "FDTrueEvRecoNumuMethod/I"
        )
        b_true_numu_nu_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoNumu", b_true_numu_nu_E, "FDTrueEvRecoNumu/F"
        )
        b_true_numu_lep_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoLepNumu", b_true_numu_lep_E, "FDTrueEvRecoLepNumu/F"
        )
        b_true_numu_had_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoHadNumu", b_true_numu_had_E, "FDTrueEvRecoHadNumu/F"
        )
        b_true_nc_nu_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoNC", b_true_nc_nu_E, "FDTrueEvRecoNC/F"
        )
        b_true_nc_lep_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoLepNC", b_true_nc_lep_E, "FDTrueEvRecoLepNC/F"
        )
        b_true_nc_had_E = array("f", [0])
        t_fdreco.Branch(
            "FDTrueEvRecoHadNC", b_true_nc_had_E, "FDTrueEvRecoHadNC/F"
        )

        # True hit info
        b_true_n_hits_Z = array("i", [0])
        t_fdreco.Branch(
            "FDTrueNumHitsZ", b_true_n_hits_Z, "FDTrueNumHitsZ/I"
        )
        b_true_sum_hits_summedadc_Z = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsSummedADCZ", b_true_sum_hits_summedadc_Z, "FDTrueSumHitsSummedADCZ/F"
        )
        b_true_sum_hits_integral_Z = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsIntegralZ", b_true_sum_hits_integral_Z, "FDTrueSumHitsIntegralZ/F"
        )
        b_true_n_hits_U = array("i", [0])
        t_fdreco.Branch(
            "FDTrueNumHitsU", b_true_n_hits_U, "FDTrueNumHitsU/I"
        )
        b_true_sum_hits_summedadc_U = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsSummedADCU", b_true_sum_hits_summedadc_U, "FDTrueSumHitsSummedADCU/F"
        )
        b_true_sum_hits_integral_U = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsIntegralU", b_true_sum_hits_integral_U, "FDTrueSumHitsIntegralU/F"
        )
        b_true_n_hits_V = array("i", [0])
        t_fdreco.Branch(
            "FDTrueNumHitsV", b_true_n_hits_V, "FDTrueNumHitsV/I"
        )
        b_true_sum_hits_summedadc_V = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsSummedADCV", b_true_sum_hits_summedadc_V, "FDTrueSumHitsSummedADCV/F"
        )
        b_true_sum_hits_integral_V = array("f", [0])
        t_fdreco.Branch(
            "FDTrueSumHitsIntegralV", b_true_sum_hits_integral_V, "FDTrueSumHitsIntegralV/F"
        )

        for i_e, (e_caf, e_eventid) in enumerate(zip(t_caf, t_eventid)):
            eventid = int(e_eventid.eventId)
            true_fd_reco = (
                f_h5["insidendlar_fd_reco"][f_h5["insidendlar_fd_reco"]["eventID"] == eventid]
            )
            pred_fd_reco = (
                f_h5["predresp_fd_reco"][f_h5["predresp_fd_reco"]["eventID"] == eventid]
            )
            b_eventid = eventid

            if len(true_fd_reco) and len(pred_fd_reco) :
                b_pred_cvn_anu[0] = float(pred_fd_reco["antinu_score"])
                b_pred_cvn_nue[0] = float(pred_fd_reco["nue_score"])
                b_pred_cvn_numu[0] = float(pred_fd_reco["numu_score"])
                b_pred_cvn_nc[0] = float(pred_fd_reco["nc_score"])
                b_pred_cvn_nutau[0] = float(pred_fd_reco["nutau_score"])
                b_pred_cvn_0p[0] = float(pred_fd_reco["0_p_score"])
                b_pred_cvn_1p[0] = float(pred_fd_reco["1_p_score"])
                b_pred_cvn_2p[0] = float(pred_fd_reco["2_p_score"])
                b_pred_cvn_np[0] = float(pred_fd_reco["N_p_score"])
                b_pred_cvn_0pi[0] = float(pred_fd_reco["0_pi_score"])
                b_pred_cvn_1pi[0] = float(pred_fd_reco["1_pi_score"])
                b_pred_cvn_2pi[0] = float(pred_fd_reco["2_pi_score"])
                b_pred_cvn_npi[0] = float(pred_fd_reco["N_pi_score"])
                b_pred_cvn_0pi0[0] = float(pred_fd_reco["0_pi0_score"])
                b_pred_cvn_1pi0[0] = float(pred_fd_reco["1_pi0_score"])
                b_pred_cvn_2pi0[0] = float(pred_fd_reco["2_pi0_score"])
                b_pred_cvn_npi0[0] = float(pred_fd_reco["N_pi0_score"])
                b_pred_cvn_0n[0] = float(pred_fd_reco["0_n_score"])
                b_pred_cvn_1n[0] = float(pred_fd_reco["1_n_score"])
                b_pred_cvn_2n[0] = float(pred_fd_reco["2_n_score"])
                b_pred_cvn_nn[0] = float(pred_fd_reco["N_n_score"])

                b_pred_numu_reco_method[0] = int(pred_fd_reco["numu_reco_method"])
                b_pred_numu_nu_E[0] = float(pred_fd_reco["numu_nu_E"])
                b_pred_numu_lep_E[0] = float(pred_fd_reco["numu_lep_E"])
                b_pred_numu_had_E[0] = float(pred_fd_reco["numu_had_E"])
                b_pred_nc_nu_E[0] = float(pred_fd_reco["nc_nu_E"])
                b_pred_nc_lep_E[0] = float(pred_fd_reco["nc_lep_E"])
                b_pred_nc_had_E[0] = float(pred_fd_reco["nc_had_E"])

                b_pred_n_hits_Z[0] = int(pred_fd_reco["n_hits_z"])
                b_pred_sum_hits_summedadc_Z[0] = float(pred_fd_reco["sum_hits_summedadc_z"])
                b_pred_sum_hits_integral_Z[0] = float(pred_fd_reco["sum_hits_integral_z"])
                b_pred_n_hits_U[0] = int(pred_fd_reco["n_hits_u"])
                b_pred_sum_hits_summedadc_U[0] = float(pred_fd_reco["sum_hits_summedadc_u"])
                b_pred_sum_hits_integral_U[0] = float(pred_fd_reco["sum_hits_integral_u"])
                b_pred_n_hits_V[0] = int(pred_fd_reco["n_hits_v"])
                b_pred_sum_hits_summedadc_V[0] = float(pred_fd_reco["sum_hits_summedadc_v"])
                b_pred_sum_hits_integral_V[0] = float(pred_fd_reco["sum_hits_integral_v"])

                b_true_cvn_anu[0] = float(true_fd_reco["antinu_score"])
                b_true_cvn_nue[0] = float(true_fd_reco["nue_score"])
                b_true_cvn_numu[0] = float(true_fd_reco["numu_score"])
                b_true_cvn_nc[0] = float(true_fd_reco["nc_score"])
                b_true_cvn_nutau[0] = float(true_fd_reco["nutau_score"])
                b_true_cvn_0p[0] = float(true_fd_reco["0_p_score"])
                b_true_cvn_1p[0] = float(true_fd_reco["1_p_score"])
                b_true_cvn_2p[0] = float(true_fd_reco["2_p_score"])
                b_true_cvn_np[0] = float(true_fd_reco["N_p_score"])
                b_true_cvn_0pi[0] = float(true_fd_reco["0_pi_score"])
                b_true_cvn_1pi[0] = float(true_fd_reco["1_pi_score"])
                b_true_cvn_2pi[0] = float(true_fd_reco["2_pi_score"])
                b_true_cvn_npi[0] = float(true_fd_reco["N_pi_score"])
                b_true_cvn_0pi0[0] = float(true_fd_reco["0_pi0_score"])
                b_true_cvn_1pi0[0] = float(true_fd_reco["1_pi0_score"])
                b_true_cvn_2pi0[0] = float(true_fd_reco["2_pi0_score"])
                b_true_cvn_npi0[0] = float(true_fd_reco["N_pi0_score"])
                b_true_cvn_0n[0] = float(true_fd_reco["0_n_score"])
                b_true_cvn_1n[0] = float(true_fd_reco["1_n_score"])
                b_true_cvn_2n[0] = float(true_fd_reco["2_n_score"])
                b_true_cvn_nn[0] = float(true_fd_reco["N_n_score"])

                b_true_numu_reco_method[0] = int(true_fd_reco["numu_reco_method"])
                b_true_numu_nu_E[0] = float(true_fd_reco["numu_nu_E"])
                b_true_numu_lep_E[0] = float(true_fd_reco["numu_lep_E"])
                b_true_numu_had_E[0] = float(true_fd_reco["numu_had_E"])
                b_true_nc_nu_E[0] = float(true_fd_reco["nc_nu_E"])
                b_true_nc_lep_E[0] = float(true_fd_reco["nc_lep_E"])
                b_true_nc_had_E[0] = float(true_fd_reco["nc_had_E"])

                b_true_n_hits_Z[0] = int(true_fd_reco["n_hits_z"])
                b_true_sum_hits_summedadc_Z[0] = float(true_fd_reco["sum_hits_summedadc_z"])
                b_true_sum_hits_integral_Z[0] = float(true_fd_reco["sum_hits_integral_z"])
                b_true_n_hits_U[0] = int(true_fd_reco["n_hits_u"])
                b_true_sum_hits_summedadc_U[0] = float(true_fd_reco["sum_hits_summedadc_u"])
                b_true_sum_hits_integral_U[0] = float(true_fd_reco["sum_hits_integral_u"])
                b_true_n_hits_V[0] = int(true_fd_reco["n_hits_v"])
                b_true_sum_hits_summedadc_V[0] = float(true_fd_reco["sum_hits_summedadc_v"])
                b_true_sum_hits_integral_V[0] = float(true_fd_reco["sum_hits_integral_v"])

            else:
                b_pred_cvn_anu[0] = -999.0
                b_pred_cvn_nue[0] = -999.0
                b_pred_cvn_numu[0] = -999.0
                b_pred_cvn_nc[0] = -999.0
                b_pred_cvn_nutau[0] = -999.0
                b_pred_cvn_0p[0] = -999.0
                b_pred_cvn_1p[0] = -999.0
                b_pred_cvn_2p[0] = -999.0
                b_pred_cvn_np[0] = -999.0
                b_pred_cvn_0pi[0] = -999.0
                b_pred_cvn_1pi[0] = -999.0
                b_pred_cvn_2pi[0] = -999.0
                b_pred_cvn_npi[0] = -999.0
                b_pred_cvn_0pi0[0] = -999.0
                b_pred_cvn_1pi0[0] = -999.0
                b_pred_cvn_2pi0[0] = -999.0
                b_pred_cvn_npi0[0] = -999.0
                b_pred_cvn_0n[0] = -999.0
                b_pred_cvn_1n[0] = -999.0
                b_pred_cvn_2n[0] = -999.0
                b_pred_cvn_nn[0] = -999.0

                b_pred_numu_reco_method[0] = -999
                b_pred_numu_nu_E[0] = -999.0
                b_pred_numu_lep_E[0] = -999.0
                b_pred_numu_had_E[0] = -999.0
                b_pred_nc_nu_E[0] = -999.0
                b_pred_nc_lep_E[0] = -999.0
                b_pred_nc_had_E[0] = -999.0

                b_pred_n_hits_Z[0] = -999
                b_pred_sum_hits_summedadc_Z[0] = -999.0
                b_pred_sum_hits_integral_Z[0] = -999.0
                b_pred_n_hits_U[0] = -999
                b_pred_sum_hits_summedadc_U[0] = -999.0
                b_pred_sum_hits_integral_U[0] = -999.0
                b_pred_n_hits_V[0] = -999
                b_pred_sum_hits_summedadc_V[0] = -999.0
                b_pred_sum_hits_integral_V[0] = -999.0

                b_true_cvn_anu[0] = -999.0
                b_true_cvn_nue[0] = -999.0
                b_true_cvn_numu[0] = -999.0
                b_true_cvn_nc[0] = -999.0
                b_true_cvn_nutau[0] = -999.0
                b_true_cvn_0p[0] = -999.0
                b_true_cvn_1p[0] = -999.0
                b_true_cvn_2p[0] = -999.0
                b_true_cvn_np[0] = -999.0
                b_true_cvn_0pi[0] = -999.0
                b_true_cvn_1pi[0] = -999.0
                b_true_cvn_2pi[0] = -999.0
                b_true_cvn_npi[0] = -999.0
                b_true_cvn_0pi0[0] = -999.0
                b_true_cvn_1pi0[0] = -999.0
                b_true_cvn_2pi0[0] = -999.0
                b_true_cvn_npi0[0] = -999.0
                b_true_cvn_0n[0] = -999.0
                b_true_cvn_1n[0] = -999.0
                b_true_cvn_2n[0] = -999.0
                b_true_cvn_nn[0] = -999.0

                b_true_numu_reco_method[0] = -999
                b_true_numu_nu_E[0] = -999.0
                b_true_numu_lep_E[0] = -999.0
                b_true_numu_had_E[0] = -999.0
                b_true_nc_nu_E[0] = -999.0
                b_true_nc_lep_E[0] = -999.0
                b_true_nc_had_E[0] = -999.0

                b_true_n_hits_Z[0] = -999
                b_true_sum_hits_summedadc_Z[0] = -999.0
                b_true_sum_hits_integral_Z[0] = -999.0
                b_true_n_hits_U[0] = -999
                b_true_sum_hits_summedadc_U[0] = -999.0
                b_true_sum_hits_integral_U[0] = -999.0
                b_true_n_hits_V[0] = -999
                b_true_sum_hits_summedadc_V[0] = -999.0
                b_true_sum_hits_integral_V[0] = -999.0

            t_fdreco.Fill()

        t_caf.AddFriend("FDRecoFriend")
        f_caf.Write()
        f_caf.Close()
        f_h5.close()

        gather_outputs(f_caf_name, f_h5_name, out_dir=args.output_dir)

""" Helpers """

def gather_inputs(f_h5_path, caf_dir):
    f_h5_name = os.path.basename(f_h5_path)
    proc = subprocess.Popen(["ifdh", "cp", f_h5_path, f_h5_name], stdout=subprocess.PIPE)
    proc.wait()

    f_prefix = ".".join(f_h5_name.split(".")[:2])
    f_caf_name = f_prefix + ".nd.CAF.root"
    f_caf_path = os.path.join(caf_dir, f_caf_name)
    proc = subprocess.Popen(["ifdh", "cp", f_caf_path, f_caf_name], stdout=subprocess.PIPE)
    proc.wait()

    return f_h5_name, f_caf_name

def gather_outputs(f_caf_name, f_h5_name, out_dir=None):
    if out_dir is not None:
        f_caf_name_out = f_caf_name.rstrip(".root") + ".fdreco.fdresppredreco.root"
        f_caf_path_out = os.path.join(out_dir, f_caf_name_out)
        proc = subprocess.Popen(["ifdh", "cp", f_caf_name, f_caf_path_out], stdout=subprocess.PIPE)
        proc.wait()

    proc = subprocess.Popen(["rm", "-v", f_h5_name], stdout=subprocess.PIPE)
    proc.wait()
    proc = subprocess.Popen(["rm", "-v", f_caf_name], stdout=subprocess.PIPE)
    proc.wait()

""" End helpers """

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("h5_dir", type=str, help="input dir with ndfd h5 files")
    parser.add_argument("caf_dir", type=str, help="input dir with nd caf files")
    parser.add_argument("output_dir", type=str, help="output dir")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main(parse_arguments())
