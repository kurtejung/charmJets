#!/bin/bash

# SSV + JP
root -l -q bfractionVsJetPtPP_CJet_v3.C\+\(\"discr_ssvHighEff\",2,0,\"SSVHE\",0,100,-2,2,0,1,0,0,1,0,1,0\)

# dR + JP
root -l -q bfractionVsJetPtPP_CJet_v3.C\+\(\"djetR\",0.3,0,\"djetR\",0,100,-2,2,0,1,0,0,0,0,1,1\)

# dR + SSV + CSV + JP
root -l -q bfractionVsJetPtPP_CJet_v3.C\+\(\"discr_ssvHighEff\",2.0,0,\"SSVHE\",0,100,-2,2,0,1,0,0,1,1,1,1\)

#do dR + SSV + JP
root -l -q bfractionVsJetPtPP_CJet_v3.C\+\(\"discr_ssvHighEff\",2.0,0,\"SSVHE\",0,100,-2,2,0,1,0,0,1,0,1,1\)
