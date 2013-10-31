#!/bin/bash


#change for signal
#declare -a sigfiles=( "radion_300" "radion_500" "radion_700" "radion_1000") 
# n files of 500 events
#declare -a MR=(300 500 700) 
#change for radion mass cut
#sed -i -r -e "s/RadMass > [0-9]+\-50/RadMass \> ${MR[$j]}-50/g" ggAnalysis_full_fit_tree.C
#sed -i -r -e "s/RadMass < [0-9]+\+50/RadMass \< ${MR[$j]}+50/g" ggAnalysis_full_fit_tree.C
#change for data file
#declare -a radion=("300" "500" "700" "1000")
#declare -a winl=("270" "470" "670" "970")
#declare -a winu=("330" "530" "730" "1030")

declare -a radion=("300" "500" "700" "1000")
declare -a winl=("285" "450" "650" "950")
declare -a winu=("315" "550" "750" "1050")

# doWP.cc
# models_mtot_exp.rs

for (( i = 4 ; i < 5 ; i++ )); do # for each working point
  mkdir radlim_CSV_WP$i
  #  for (( j = 0 ; j < ${#radion[@]} ; j++ )); do # for each mass 
  for (( j = 1 ; j < 4 ; j++ )); do # for each mass 
  #  for (( j = 0 ; j < 1 ; j++ )); do # for each mass 
	# create the datacard and the workspaces
	#check the name on R2GGBBFitter.cc!!   legmc->SetHeader("300 GeV | CIC + X jets selection");
        # TString ssignal   = "./MiniTrees/jetwin/WP1/1000.root";
        # TString ddata   = "./MiniTrees/jetwin/WP1/Data_1000.root";
	# things to choose on the runcard
	# the signal file
	sed -i -r -e "s/WP[0-9]/WP$i/g" doWP.cc
	sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" doWP.cc
	sed -i -r -e "s/m[0-9]+_8TeV_nm_m[0-9]+.root/m${radion[$j]}_8TeV_nm_m${radion[$j]}.root/g" doWP.cc

	#sed -i -r -e "s/[0-9]+\_minimal.root/${radion[$j]}\_minimal.root/g" doWP.cc
	#sed -i -r -e "s/[0-9]+\_regression/${radion[$j]}\_regression/g" doWP.cc
#	sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" doWP.cc
#	sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" models_mtot_exp.rs
	# the legend
	sed -i -r -e "s/[0-9]+ GeV\"\);/${radion[$j]} GeV\"\);/g" doWP.cc
        # name of files
	sed -i -r -e "s/mH[0-9]+/mH${radion[$j]}/g" doWP.cc
	# the window mtot > 550 || mtot < 450
	#sed -i -r -e "s/mtot > [0-9]+ \|\| mtot < [0-9]+/mtot > ${winu[$j]} \|\| mtot < ${winl[$j]}/g" doWP.cc
echo WP$i MR ${radion[$j]}
	# the mass to fit around
	sed -i -r -e "s/runfits\([0-9]+/runfits\(${radion[$j]}/g" runfits.C
	# the fit model
	# mtot_sig_m0[500.0, 450, 550];
	# mtot_sig_m0_cat0[500.0, 450, 550];
	# mtot_sig_m0_cat1[500.0, 450, 550];
	# mtot_sig_m0_cat2[500.0, 450, 550];
	sed -i -r -e "s/mtot_sig_m0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	sed -i -r -e "s/mtot_sig_m0_cat0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	sed -i -r -e "s/mtot_sig_m0_cat1\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat1\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	#sed -i -r -e "s/mtot_sig_m0_cat2\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat2\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_exp.rs
	mkdir radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	root -l -q runfits.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}.txt
	cp hgg.* radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	cp .txt radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	#also colect the plots
	mv databkgoversig*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv sigmodel*  radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	mv remenber.txt radlim${radion[$j]}_CSV
	## create limits root files for each mass
#	cd radlim${radion[$j]}_CSV
	cd radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV_onecat.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecat.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecat.root
	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
#	combine -M Asymptotic hgg.mH${radion[$j]}_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.txt
	#combine hgg.mH${radion[$j]}.0_8TeV.txt -M Asymptotic -m ${radion[$j]} -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
        #combine hgg.mH${radion[$j]}.0_8TeVonecat.txt -M Asymptotic -m 125 -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}onecat_nosyst.txt
	#mv higgsCombineTest.Asymptotic.mH125.root higgsCombineTest.Asymptotic.mR${radion[$j]}_onecat_nosyst.root #adapt this output name to the cut
	#mv higgsCombineTest.Asymptotic.mH125.root higgsCombineTest.Asymptotic.mR${radion[$j]}_CSV.root #adapt this output name to the cut
	cd ../..
  done # mass
#  sed -i -r -e "s/WP[0-9]/WP$i/g" Brazilianflag.cc
#  root -l -q Brazilianflag.cc
done # wp
#sed -i -r -e "s/mR300.Xcut/Mcut_cutbased/g" Brazilianflag.cc #adapt this name to the cut
#root -l Brazilianflag.cc

