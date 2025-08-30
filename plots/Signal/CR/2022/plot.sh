#!/bin/bash
python CR_2022.py --hist-name DileptonMass_Central --output-name DileptonMass_CR_2022_EE --x-range 50,500  --logy --ymax 1e3 --ymin 1e-3 --x-title "m_{ll} [GeV]" --custom-bins "50,70,90,130,500"
python CR_2022.py --hist-name LeadingBJetPt_Central --output-name LeadingBJetpt_CR_2022_EE --x-range 30,500  --logy --ymax 1e3 --ymin 1e-3 --x-title "B-jet p_{T} [GeV]" --custom-bins "30,40,50,70,90,120,150,200,300,5000"
python CR_2022.py --hist-name LeadingMuonPt_Central --output-name LeadingMuonpt_CR_2022_EE --x-range 50,500  --logy --ymax 1e3 --ymin 1e-3 --x-title "Leading Muon p_{T} [GeV]" --custom-bins "50,70,90,110,130,150,200,5000"
python CR_2022.py --hist-name LeadingTopJetPt_Central --output-name LeadingTopJetpt_CR_2022_EE --x-range 300,1000  --logy --ymax 1e3 --ymin 1e-3 --x-title "Top Jet p_{T} [GeV]" --custom-bins "300,350,400,450,500,550,600,650,5000"
python CR_2022.py --hist-name SubleadingMuonPt_Central --output-name SubleadingMuonpt_CR_2022_EE --x-range 0,200  --logy --ymax 1e3 --ymin 1e-3 --x-title "Subleading Muon p_{T} [GeV]" --custom-bins "0,10,20,50,100,5000"
python CR_2022.py --hist-name Topjetnum --output-name TopJetnum_CR_2022_EE --x-range 0,5  --logy --ymax 1e4 --ymin 1e-3 --x-title "Topjetnum"
python CR_2022.py --hist-name WRMass_Central --output-name WRMass_CR_2022_EE --x-range 200,2000  --logy --ymax 1e3 --ymin 1e-3 --x-title "m_{WR} [GeV]" --custom-bins "200,300,400,500,600,700,800,900,1000,1200,1500,2000,3000,5000"
python CR_2022.py --hist-name Bjetnum --output-name Bjetnum_CR_2022_EE --x-range 0,10  --logy --ymax 1e4 --ymin 1e-3 --x-title "Bjetnum"
