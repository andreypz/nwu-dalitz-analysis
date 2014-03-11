#ifndef PHOTONSCALECORRECTIONS_HH
#define PHOTONSCALECORRECTIONS_HH

#include <iostream>
#include <cassert>
#include <TRandom3.h>

Double_t correctedPhotonEnergy( Double_t inputEnergy, Double_t scEta, Double_t R9, Int_t run, Int_t type, string targetDataset, Bool_t isMC, TRandom3 *random ) {

  double scale = 1.0;
  //********************************************************
  //Std Ecal Energy
  //********************************************************
  if (type == 0) {

    if (targetDataset == "HCP2012") {

      if ( fabs(scEta) < 1 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 1.0057;
        if (run>=190782 && run<=191042) scale = 1.0115;
        if (run>=191043 && run<=193555) scale = 1.0029;
        if (run>=193556 && run<=194150) scale = 1.0018;
        if (run>=194151 && run<=194532) scale = 1.0016;
        if (run>=194533 && run<=195113) scale = 1.0015;
        if (run>=195114 && run<=195915) scale = 1.0015;
        if (run>=195916 && run<=198115) scale = 1.0006;
        if (run>=198116 && run<=199803) scale = 1.0046;
        if (run>=199804 && run<=200048) scale = 1.0053;
        if (run>=200049 && run<=200151) scale = 1.0064;
        if (run>=200152 && run<=200490) scale = 1.0045;
        if (run>=200491 && run<=200531) scale = 1.0058;
        if (run>=200532 && run<=201656) scale = 1.0045;
        if (run>=201657 && run<=202305) scale = 1.0053;
        if (run>=202305 && run<=203002) scale = 1.0065;
      } else if (fabs(scEta) < 1 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 1.0042;
        if (run>=190782 && run<=191042) scale = 1.0099;
        if (run>=191043 && run<=193555) scale = 1.0014;
        if (run>=193556 && run<=194150) scale = 1.0002;
        if (run>=194151 && run<=194532) scale = 1.0001;
        if (run>=194533 && run<=195113) scale = 1.0000;
        if (run>=195114 && run<=195915) scale = 0.9999;
        if (run>=195916 && run<=198115) scale = 0.9991;
        if (run>=198116 && run<=199803) scale = 1.0031;
        if (run>=199804 && run<=200048) scale = 1.0037;
        if (run>=200049 && run<=200151) scale = 1.0049;
        if (run>=200152 && run<=200490) scale = 1.0029;
        if (run>=200491 && run<=200531) scale = 1.0042;
        if (run>=200532 && run<=201656) scale = 1.0030;
        if (run>=201657 && run<=202305) scale = 1.0038;
        if (run>=202305 && run<=203002) scale = 1.0050;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 1.0048;
        if (run>=190782 && run<=191042) scale = 1.0070;
        if (run>=191043 && run<=193555) scale = 1.0017;
        if (run>=193556 && run<=194150) scale = 0.9979;
        if (run>=194151 && run<=194532) scale = 0.9980;
        if (run>=194533 && run<=195113) scale = 0.9993;
        if (run>=195114 && run<=195915) scale = 0.9981;
        if (run>=195916 && run<=198115) scale = 0.9971;
        if (run>=198116 && run<=199803) scale = 1.0014;
        if (run>=199804 && run<=200048) scale = 1.0020;
        if (run>=200049 && run<=200151) scale = 1.0026;
        if (run>=200152 && run<=200490) scale = 1.0012;
        if (run>=200491 && run<=200531) scale = 1.0014;
        if (run>=200532 && run<=201656) scale = 1.0002;
        if (run>=201657 && run<=202305) scale = 1.0002;
        if (run>=202305 && run<=203002) scale = 1.0010;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9953;
        if (run>=190782 && run<=191042) scale = 0.9975;
        if (run>=191043 && run<=193555) scale = 0.9921;
        if (run>=193556 && run<=194150) scale = 0.9883;
        if (run>=194151 && run<=194532) scale = 0.9884;
        if (run>=194533 && run<=195113) scale = 0.9897;
        if (run>=195114 && run<=195915) scale = 0.9884;
        if (run>=195916 && run<=198115) scale = 0.9875;
        if (run>=198116 && run<=199803) scale = 0.9918;
        if (run>=199804 && run<=200048) scale = 0.9924;
        if (run>=200049 && run<=200151) scale = 0.9930;
        if (run>=200152 && run<=200490) scale = 0.9916;
        if (run>=200491 && run<=200531) scale = 0.9918;
        if (run>=200532 && run<=201656) scale = 0.9906;
        if (run>=201657 && run<=202305) scale = 0.9906;
        if (run>=202305 && run<=203002) scale = 0.9914;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 0.9953;
        if (run>=190782 && run<=191042) scale = 1.0006;
        if (run>=191043 && run<=193555) scale = 0.9999;
        if (run>=193556 && run<=194150) scale = 1.0009;
        if (run>=194151 && run<=194532) scale = 1.0019;
        if (run>=194533 && run<=195113) scale = 1.0018;
        if (run>=195114 && run<=195915) scale = 1.0019;
        if (run>=195916 && run<=198115) scale = 1.0035;
        if (run>=198116 && run<=199803) scale = 0.9982;
        if (run>=199804 && run<=200048) scale = 0.9972;
        if (run>=200049 && run<=200151) scale = 1.0001;
        if (run>=200152 && run<=200490) scale = 1.0002;
        if (run>=200491 && run<=200531) scale = 0.9971;
        if (run>=200532 && run<=201656) scale = 0.9978;
        if (run>=201657 && run<=202305) scale = 0.9984;
        if (run>=202305 && run<=203002) scale = 1.0003;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9878;
        if (run>=190782 && run<=191042) scale = 0.9931;
        if (run>=191043 && run<=193555) scale = 0.9924;
        if (run>=193556 && run<=194150) scale = 0.9934;
        if (run>=194151 && run<=194532) scale = 0.9945;
        if (run>=194533 && run<=195113) scale = 0.9943;
        if (run>=195114 && run<=195915) scale = 0.9945;
        if (run>=195916 && run<=198115) scale = 0.9961;
        if (run>=198116 && run<=199803) scale = 0.9907;
        if (run>=199804 && run<=200048) scale = 0.9898;
        if (run>=200049 && run<=200151) scale = 0.9927;
        if (run>=200152 && run<=200490) scale = 0.9928;
        if (run>=200491 && run<=200531) scale = 0.9896;
        if (run>=200532 && run<=201656) scale = 0.9903;
        if (run>=201657 && run<=202305) scale = 0.9910;
        if (run>=202305 && run<=203002) scale = 0.9928;
      } else if (fabs(scEta) >=2 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 1.0032;
        if (run>=190782 && run<=191042) scale = 1.0085;
        if (run>=191043 && run<=193555) scale = 1.0074;
        if (run>=193556 && run<=194150) scale = 1.0082;
        if (run>=194151 && run<=194532) scale = 1.0086;
        if (run>=194533 && run<=195113) scale = 1.0069;
        if (run>=195114 && run<=195915) scale = 1.0070;
        if (run>=195916 && run<=198115) scale = 1.0066;
        if (run>=198116 && run<=199803) scale = 1.0127;
        if (run>=199804 && run<=200048) scale = 1.0132;
        if (run>=200049 && run<=200151) scale = 1.0161;
        if (run>=200152 && run<=200490) scale = 1.0155;
        if (run>=200491 && run<=200531) scale = 1.0161;
        if (run>=200532 && run<=201656) scale = 1.0138;
        if (run>=201657 && run<=202305) scale = 1.0182;
        if (run>=202305 && run<=203002) scale = 1.0210;
      } else if (fabs(scEta) >= 2 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9897;
        if (run>=190782 && run<=191042) scale = 0.9951;
        if (run>=191043 && run<=193555) scale = 0.9939;
        if (run>=193556 && run<=194150) scale = 0.9947;
        if (run>=194151 && run<=194532) scale = 0.9952;
        if (run>=194533 && run<=195113) scale = 0.9935;
        if (run>=195114 && run<=195915) scale = 0.9935;
        if (run>=195916 && run<=198115) scale = 0.9931;
        if (run>=198116 && run<=199803) scale = 0.9993;
        if (run>=199804 && run<=200048) scale = 0.9998;
        if (run>=200049 && run<=200151) scale = 1.0028;
        if (run>=200152 && run<=200490) scale = 1.0021;
        if (run>=200491 && run<=200531) scale = 1.0028;
        if (run>=200532 && run<=201656) scale = 1.0004;
        if (run>=201657 && run<=202305) scale = 1.0049;
        if (run>=202305 && run<=203002) scale = 1.0077;
      }
    } //end if dataset == HCP2012 ( j13,a06,a24 rereco + pr-v2 )

    else if (targetDataset == "ICHEP2012") {
      // values from https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/873.html	
      if ( fabs(scEta) < 1 && R9<0.94) {
        if (run>=190450 && run<=190781) scale = 1.0021; 
        if (run>=190782 && run<=190949) scale = 1.0154;   
        if (run>=190950 && run<=191833) scale = 1.0046;   
        if (run>=191834 && run<=193686) scale = 1.0017;   
        if (run>=193746 && run<=194210) scale = 1.0020;   
        if (run>=194211 && run<=194479) scale = 1.0037;   
        if (run>=194480 && run<=195147) scale = 1.0047;   
        if (run>=195148 && run<=195350) scale = 1.0053;   
        if (run>=195396 && run<=195530) scale = 1.0042;   
        if (run>=195531 && run<=196531) scale = 0.9981;   
      } else if (fabs(scEta) < 1 && R9>=0.94) {
        if (run>=190450 && run<=190781) scale = 0.9962; 
        if (run>=190782 && run<=190949) scale = 1.0096;   
        if (run>=190950 && run<=191833) scale = 0.9988;   
        if (run>=191834 && run<=193686) scale = 0.9958;   
        if (run>=193746 && run<=194210) scale = 0.9962;   
        if (run>=194211 && run<=194479) scale = 0.9979;   
        if (run>=194480 && run<=195147) scale = 0.9989;   
        if (run>=195148 && run<=195350) scale = 0.9995;   
        if (run>=195396 && run<=195530) scale = 0.9984;   
        if (run>=195531 && run<=196531) scale = 0.9922;   
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
        if (run>=190450 && run<=190781) scale = 1.0133; 
        if (run>=190782 && run<=190949) scale = 0.9997;   
        if (run>=190950 && run<=191833) scale = 1.0134;   
        if (run>=191834 && run<=193686) scale = 1.0104;   
        if (run>=193746 && run<=194210) scale = 1.0094;   
        if (run>=194211 && run<=194479) scale = 1.0118;   
        if (run>=194480 && run<=195147) scale = 1.0137;   
        if (run>=195148 && run<=195350) scale = 1.0142;   
        if (run>=195396 && run<=195530) scale = 1.0129;   
        if (run>=195531 && run<=196531) scale = 1.0065;   
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
        if (run>=190450 && run<=190781) scale = 1.0020; 
        if (run>=190782 && run<=190949) scale = 0.9883;   
        if (run>=190950 && run<=191833) scale = 1.0021;   
        if (run>=191834 && run<=193686) scale = 0.9991;   
        if (run>=193746 && run<=194210) scale = 0.9980;   
        if (run>=194211 && run<=194479) scale = 1.0005;   
        if (run>=194480 && run<=195147) scale = 1.0024;   
        if (run>=195148 && run<=195350) scale = 1.0029;   
        if (run>=195396 && run<=195530) scale = 1.0016;   
        if (run>=195531 && run<=196531) scale = 0.9951;   
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
        if (run>=190450 && run<=190781) scale = 0.9989; 
        if (run>=190782 && run<=190949) scale = 1.0123;   
        if (run>=190950 && run<=191833) scale = 1.0042;   
        if (run>=191834 && run<=193686) scale = 1.0037;   
        if (run>=193746 && run<=194210) scale = 1.0047;   
        if (run>=194211 && run<=194479) scale = 1.0037;   
        if (run>=194480 && run<=195147) scale = 1.0030;   
        if (run>=195148 && run<=195350) scale = 1.0051;   
        if (run>=195396 && run<=195530) scale = 1.0003;   
        if (run>=195531 && run<=196531) scale = 1.0052;   
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
        if (run>=190450 && run<=190781) scale = 0.9931; 
        if (run>=190782 && run<=190949) scale = 1.0066;   
        if (run>=190950 && run<=191833) scale = 0.9985;   
        if (run>=191834 && run<=193686) scale = 0.9979;   
        if (run>=193746 && run<=194210) scale = 0.9990;   
        if (run>=194211 && run<=194479) scale = 0.9979;   
        if (run>=194480 && run<=195147) scale = 0.9972;   
        if (run>=195148 && run<=195350) scale = 0.9994;   
        if (run>=195396 && run<=195530) scale = 0.9945;   
        if (run>=195531 && run<=196531) scale = 0.9994;   
      } else if (fabs(scEta) >=2 && R9<0.94) {
        if (run>=190450 && run<=190781) scale = 1.0052; 
        if (run>=190782 && run<=190949) scale = 1.0077;   
        if (run>=190950 && run<=191833) scale = 0.9900;   
        if (run>=191834 && run<=193686) scale = 0.9893;   
        if (run>=193746 && run<=194210) scale = 1.0042;   
        if (run>=194211 && run<=194479) scale = 1.0036;   
        if (run>=194480 && run<=195147) scale = 1.0069;   
        if (run>=195148 && run<=195350) scale = 1.0133;   
        if (run>=195396 && run<=195530) scale = 0.9999;   
        if (run>=195531 && run<=196531) scale = 1.0475;   
      } else if (fabs(scEta) >= 2 && R9>=0.94) {
        if (run>=190450 && run<=190781) scale = 0.9899; 
        if (run>=190782 && run<=190949) scale = 0.9924;   
        if (run>=190950 && run<=191833) scale = 0.9745;   
        if (run>=191834 && run<=193686) scale = 0.9738;   
        if (run>=193746 && run<=194210) scale = 0.9889;   
        if (run>=194211 && run<=194479) scale = 0.9883;   
        if (run>=194480 && run<=195147) scale = 0.9916;   
        if (run>=195148 && run<=195350) scale = 0.9982;   
        if (run>=195396 && run<=195530) scale = 0.9845;   
        if (run>=195531 && run<=196531) scale = 1.0329;   
      }
    } //end if dataset == ICHEP2012 ( prompt reco )

    if (targetDataset == "Moriond2014") {
      if ( fabs(scEta) < 1 && R9<0.94) {
	if (run>=190645 && run<=190781) scale = 0.9999;
	if (run>=190782 && run<=191042) scale = 1.0064;
	if (run>=191043 && run<=191720) scale = 1.0009;
	if (run>=191721 && run<=193833) scale = 1.0001;
	if (run>=193834 && run<=194116) scale = 1.0008;
	if (run>=194117 && run<=194427) scale = 1.0008;
	if (run>=194428 && run<=194618) scale = 1.0004;
	if (run>=194619 && run<=194789) scale = 1.0008;
	if (run>=194790 && run<=195111) scale = 1.0006;
	if (run>=195112 && run<=195377) scale = 1.0010;
	if (run>=195378 && run<=195398) scale = 1.0007;
	if (run>=195399 && run<=195657) scale = 1.0009;
	if (run>=195658 && run<=195918) scale = 1.0011;
	if (run>=195919 && run<=196198) scale = 1.0008;
	if (run>=196199 && run<=196356) scale = 1.0011;
	if (run>=196357 && run<=198115) scale = 1.0010;
	if (run>=198116 && run<=198940) scale = 1.0007;
	if (run>=198941 && run<=199317) scale = 1.0009;
	if (run>=199318 && run<=199428) scale = 1.0004;
	if (run>=199429 && run<=199697) scale = 1.0012;
	if (run>=199698 && run<=199832) scale = 1.0009;
	if (run>=199833 && run<=199960) scale = 1.0011;
	if (run>=199961 && run<=200151) scale = 1.0012;
	if (run>=200152 && run<=200490) scale = 1.0010;
	if (run>=200491 && run<=200991) scale = 1.0014;
	if (run>=200992 && run<=201201) scale = 1.0007;
	if (run>=201202 && run<=201624) scale = 1.0010;
	if (run>=201625 && run<=201707) scale = 1.0012;
	if (run>=201708 && run<=202059) scale = 1.0013;
	if (run>=202060 && run<=202204) scale = 1.0016;
	if (run>=202205 && run<=202332) scale = 1.0015;
	if (run>=202333 && run<=202972) scale = 1.0016;
	if (run>=202973 && run<=203002) scale = 1.0016;
	if (run>=203003 && run<=203852) scale = 1.0028;
	if (run>=203853 && run<=204099) scale = 1.0004;
	if (run>=204100 && run<=204562) scale = 1.0002;
	if (run>=204563 && run<=205085) scale = 1.0003;
	if (run>=205086 && run<=205310) scale = 1.0002;
	if (run>=205311 && run<=205617) scale = 1.0006;
	if (run>=205618 && run<=205825) scale = 1.0006;
	if (run>=205826 && run<=206207) scale = 1.0006;
	if (run>=206208 && run<=206389) scale = 1.0011;
	if (run>=206390 && run<=206483) scale = 1.0006;
	if (run>=206484 && run<=206597) scale = 1.0004;
	if (run>=206598 && run<=206896) scale = 1.0008;
	if (run>=206897 && run<=207220) scale = 1.0014;
	if (run>=207221 && run<=207315) scale = 1.0012;
	if (run>=207316 && run<=207489) scale = 1.0009;
	if (run>=207490 && run<=207919) scale = 1.0015;
	if (run>=207920 && run<=208351) scale = 1.0010;
	if (run>=208352 && run<=208686) scale = 1.0011;
      } else if (fabs(scEta) < 1 && R9>=0.94) {
	if (run>=190645 && run<=190781) scale = 0.9977;
	if (run>=190782 && run<=191042) scale = 1.0042;
	if (run>=191043 && run<=191720) scale = 0.9987;
	if (run>=191721 && run<=193833) scale = 0.9979;
	if (run>=193834 && run<=194116) scale = 0.9986;
	if (run>=194117 && run<=194427) scale = 0.9986;
	if (run>=194428 && run<=194618) scale = 0.9982;
	if (run>=194619 && run<=194789) scale = 0.9986;
	if (run>=194790 && run<=195111) scale = 0.9984;
	if (run>=195112 && run<=195377) scale = 0.9988;
	if (run>=195378 && run<=195398) scale = 0.9985;
	if (run>=195399 && run<=195657) scale = 0.9987;
	if (run>=195658 && run<=195918) scale = 0.9989;
	if (run>=195919 && run<=196198) scale = 0.9986;
	if (run>=196199 && run<=196356) scale = 0.9989;
	if (run>=196357 && run<=198115) scale = 0.9988;
	if (run>=198116 && run<=198940) scale = 0.9985;
	if (run>=198941 && run<=199317) scale = 0.9987;
	if (run>=199318 && run<=199428) scale = 0.9982;
	if (run>=199429 && run<=199697) scale = 0.9990;
	if (run>=199698 && run<=199832) scale = 0.9987;
	if (run>=199833 && run<=199960) scale = 0.9989;
	if (run>=199961 && run<=200151) scale = 0.9990;
	if (run>=200152 && run<=200490) scale = 0.9988;
	if (run>=200491 && run<=200991) scale = 0.9992;
	if (run>=200992 && run<=201201) scale = 0.9985;
	if (run>=201202 && run<=201624) scale = 0.9988;
	if (run>=201625 && run<=201707) scale = 0.9990;
	if (run>=201708 && run<=202059) scale = 0.9991;
	if (run>=202060 && run<=202204) scale = 0.9994;
	if (run>=202205 && run<=202332) scale = 0.9993;
	if (run>=202333 && run<=202972) scale = 0.9994;
	if (run>=202973 && run<=203002) scale = 0.9994;
	if (run>=203003 && run<=203852) scale = 1.0006;
	if (run>=203853 && run<=204099) scale = 0.9982;
	if (run>=204100 && run<=204562) scale = 0.9980;
	if (run>=204563 && run<=205085) scale = 0.9981;
	if (run>=205086 && run<=205310) scale = 0.9980;
	if (run>=205311 && run<=205617) scale = 0.9984;
	if (run>=205618 && run<=205825) scale = 0.9984;
	if (run>=205826 && run<=206207) scale = 0.9984;
	if (run>=206208 && run<=206389) scale = 0.9989;
	if (run>=206390 && run<=206483) scale = 0.9984;
	if (run>=206484 && run<=206597) scale = 0.9982;
	if (run>=206598 && run<=206896) scale = 0.9986;
	if (run>=206897 && run<=207220) scale = 0.9992;
	if (run>=207221 && run<=207315) scale = 0.9990;
	if (run>=207316 && run<=207489) scale = 0.9987;
	if (run>=207490 && run<=207919) scale = 0.9993;
	if (run>=207920 && run<=208351) scale = 0.9988;
	if (run>=208352 && run<=208686) scale = 0.9989;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
	if (run>=190645 && run<=190781) scale = 1.0068;
	if (run>=190782 && run<=191042) scale = 1.0095;
	if (run>=191043 && run<=191720) scale = 1.0056;
	if (run>=191721 && run<=193833) scale = 1.0080;
	if (run>=193834 && run<=194116) scale = 1.0071;
	if (run>=194117 && run<=194427) scale = 1.0064;
	if (run>=194428 && run<=194618) scale = 1.0072;
	if (run>=194619 && run<=194789) scale = 1.0073;
	if (run>=194790 && run<=195111) scale = 1.0078;
	if (run>=195112 && run<=195377) scale = 1.0067;
	if (run>=195378 && run<=195398) scale = 1.0065;
	if (run>=195399 && run<=195657) scale = 1.0082;
	if (run>=195658 && run<=195918) scale = 1.0073;
	if (run>=195919 && run<=196198) scale = 1.0069;
	if (run>=196199 && run<=196356) scale = 1.0073;
	if (run>=196357 && run<=198115) scale = 1.0070;
	if (run>=198116 && run<=198940) scale = 1.0072;
	if (run>=198941 && run<=199317) scale = 1.0067;
	if (run>=199318 && run<=199428) scale = 1.0070;
	if (run>=199429 && run<=199697) scale = 1.0078;
	if (run>=199698 && run<=199832) scale = 1.0082;
	if (run>=199833 && run<=199960) scale = 1.0073;
	if (run>=199961 && run<=200151) scale = 1.0062;
	if (run>=200152 && run<=200490) scale = 1.0071;
	if (run>=200491 && run<=200991) scale = 1.0056;
	if (run>=200992 && run<=201201) scale = 1.0060;
	if (run>=201202 && run<=201624) scale = 1.0077;
	if (run>=201625 && run<=201707) scale = 1.0068;
	if (run>=201708 && run<=202059) scale = 1.0072;
	if (run>=202060 && run<=202204) scale = 1.0071;
	if (run>=202205 && run<=202332) scale = 1.0074;
	if (run>=202333 && run<=202972) scale = 1.0079;
	if (run>=202973 && run<=203002) scale = 1.0061;
	if (run>=203003 && run<=203852) scale = 1.0012;
	if (run>=203853 && run<=204099) scale = 1.0079;
	if (run>=204100 && run<=204562) scale = 1.0078;
	if (run>=204563 && run<=205085) scale = 1.0083;
	if (run>=205086 && run<=205310) scale = 1.0067;
	if (run>=205311 && run<=205617) scale = 1.0081;
	if (run>=205618 && run<=205825) scale = 1.0066;
	if (run>=205826 && run<=206207) scale = 1.0082;
	if (run>=206208 && run<=206389) scale = 1.0077;
	if (run>=206390 && run<=206483) scale = 1.0071;
	if (run>=206484 && run<=206597) scale = 1.0061;
	if (run>=206598 && run<=206896) scale = 1.0080;
	if (run>=206897 && run<=207220) scale = 1.0072;
	if (run>=207221 && run<=207315) scale = 1.0083;
	if (run>=207316 && run<=207489) scale = 1.0076;
	if (run>=207490 && run<=207919) scale = 1.0086;
	if (run>=207920 && run<=208351) scale = 1.0084;
	if (run>=208352 && run<=208686) scale = 1.0089;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
	if (run>=190645 && run<=190781) scale = 0.9947;
	if (run>=190782 && run<=191042) scale = 0.9974;
	if (run>=191043 && run<=191720) scale = 0.9935;
	if (run>=191721 && run<=193833) scale = 0.9959;
	if (run>=193834 && run<=194116) scale = 0.9950;
	if (run>=194117 && run<=194427) scale = 0.9944;
	if (run>=194428 && run<=194618) scale = 0.9951;
	if (run>=194619 && run<=194789) scale = 0.9952;
	if (run>=194790 && run<=195111) scale = 0.9957;
	if (run>=195112 && run<=195377) scale = 0.9947;
	if (run>=195378 && run<=195398) scale = 0.9944;
	if (run>=195399 && run<=195657) scale = 0.9961;
	if (run>=195658 && run<=195918) scale = 0.9953;
	if (run>=195919 && run<=196198) scale = 0.9948;
	if (run>=196199 && run<=196356) scale = 0.9952;
	if (run>=196357 && run<=198115) scale = 0.9949;
	if (run>=198116 && run<=198940) scale = 0.9951;
	if (run>=198941 && run<=199317) scale = 0.9946;
	if (run>=199318 && run<=199428) scale = 0.9949;
	if (run>=199429 && run<=199697) scale = 0.9957;
	if (run>=199698 && run<=199832) scale = 0.9961;
	if (run>=199833 && run<=199960) scale = 0.9953;
	if (run>=199961 && run<=200151) scale = 0.9941;
	if (run>=200152 && run<=200490) scale = 0.9950;
	if (run>=200491 && run<=200991) scale = 0.9935;
	if (run>=200992 && run<=201201) scale = 0.9939;
	if (run>=201202 && run<=201624) scale = 0.9956;
	if (run>=201625 && run<=201707) scale = 0.9947;
	if (run>=201708 && run<=202059) scale = 0.9951;
	if (run>=202060 && run<=202204) scale = 0.9950;
	if (run>=202205 && run<=202332) scale = 0.9953;
	if (run>=202333 && run<=202972) scale = 0.9958;
	if (run>=202973 && run<=203002) scale = 0.9940;
	if (run>=203003 && run<=203852) scale = 0.9892;
	if (run>=203853 && run<=204099) scale = 0.9958;
	if (run>=204100 && run<=204562) scale = 0.9957;
	if (run>=204563 && run<=205085) scale = 0.9962;
	if (run>=205086 && run<=205310) scale = 0.9946;
	if (run>=205311 && run<=205617) scale = 0.9960;
	if (run>=205618 && run<=205825) scale = 0.9945;
	if (run>=205826 && run<=206207) scale = 0.9961;
	if (run>=206208 && run<=206389) scale = 0.9956;
	if (run>=206390 && run<=206483) scale = 0.9950;
	if (run>=206484 && run<=206597) scale = 0.9941;
	if (run>=206598 && run<=206896) scale = 0.9959;
	if (run>=206897 && run<=207220) scale = 0.9951;
	if (run>=207221 && run<=207315) scale = 0.9962;
	if (run>=207316 && run<=207489) scale = 0.9955;
	if (run>=207490 && run<=207919) scale = 0.9965;
	if (run>=207920 && run<=208351) scale = 0.9963;
	if (run>=208352 && run<=208686) scale = 0.9968;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
	if (run>=190645 && run<=190781) scale = 1.0033;
	if (run>=190782 && run<=191042) scale = 1.0041;
	if (run>=191043 && run<=191720) scale = 1.0062;
	if (run>=191721 && run<=193833) scale = 1.0048;
	if (run>=193834 && run<=194116) scale = 1.0036;
	if (run>=194117 && run<=194427) scale = 1.0030;
	if (run>=194428 && run<=194618) scale = 1.0030;
	if (run>=194619 && run<=194789) scale = 1.0060;
	if (run>=194790 && run<=195111) scale = 1.0070;
	if (run>=195112 && run<=195377) scale = 1.0049;
	if (run>=195378 && run<=195398) scale = 1.0049;
	if (run>=195399 && run<=195657) scale = 1.0041;
	if (run>=195658 && run<=195918) scale = 1.0070;
	if (run>=195919 && run<=196198) scale = 1.0053;
	if (run>=196199 && run<=196356) scale = 1.0048;
	if (run>=196357 && run<=198115) scale = 1.0034;
	if (run>=198116 && run<=198940) scale = 1.0060;
	if (run>=198941 && run<=199317) scale = 1.0035;
	if (run>=199318 && run<=199428) scale = 1.0034;
	if (run>=199429 && run<=199697) scale = 1.0040;
	if (run>=199698 && run<=199832) scale = 1.0025;
	if (run>=199833 && run<=199960) scale = 1.0051;
	if (run>=199961 && run<=200151) scale = 1.0032;
	if (run>=200152 && run<=200490) scale = 1.0051;
	if (run>=200491 && run<=200991) scale = 1.0032;
	if (run>=200992 && run<=201201) scale = 1.0024;
	if (run>=201202 && run<=201624) scale = 1.0023;
	if (run>=201625 && run<=201707) scale = 1.0021;
	if (run>=201708 && run<=202059) scale = 1.0042;
	if (run>=202060 && run<=202204) scale = 1.0027;
	if (run>=202205 && run<=202332) scale = 1.0051;
	if (run>=202333 && run<=202972) scale = 1.0026;
	if (run>=202973 && run<=203002) scale = 1.0039;
	if (run>=203003 && run<=203852) scale = 0.9884;
	if (run>=203853 && run<=204099) scale = 1.0032;
	if (run>=204100 && run<=204562) scale = 1.0032;
	if (run>=204563 && run<=205085) scale = 1.0045;
	if (run>=205086 && run<=205310) scale = 1.0028;
	if (run>=205311 && run<=205617) scale = 1.0029;
	if (run>=205618 && run<=205825) scale = 1.0032;
	if (run>=205826 && run<=206207) scale = 1.0033;
	if (run>=206208 && run<=206389) scale = 1.0013;
	if (run>=206390 && run<=206483) scale = 1.0054;
	if (run>=206484 && run<=206597) scale = 1.0038;
	if (run>=206598 && run<=206896) scale = 1.0027;
	if (run>=206897 && run<=207220) scale = 1.0043;
	if (run>=207221 && run<=207315) scale = 1.0039;
	if (run>=207316 && run<=207489) scale = 1.0034;
	if (run>=207490 && run<=207919) scale = 1.0037;
	if (run>=207920 && run<=208351) scale = 1.0032;
	if (run>=208352 && run<=208686) scale = 1.0032;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
	if (run>=190645 && run<=190781) scale = 0.9929;
	if (run>=190782 && run<=191042) scale = 0.9937;
	if (run>=191043 && run<=191720) scale = 0.9957;
	if (run>=191721 && run<=193833) scale = 0.9943;
	if (run>=193834 && run<=194116) scale = 0.9932;
	if (run>=194117 && run<=194427) scale = 0.9926;
	if (run>=194428 && run<=194618) scale = 0.9926;
	if (run>=194619 && run<=194789) scale = 0.9955;
	if (run>=194790 && run<=195111) scale = 0.9965;
	if (run>=195112 && run<=195377) scale = 0.9945;
	if (run>=195378 && run<=195398) scale = 0.9944;
	if (run>=195399 && run<=195657) scale = 0.9937;
	if (run>=195658 && run<=195918) scale = 0.9966;
	if (run>=195919 && run<=196198) scale = 0.9948;
	if (run>=196199 && run<=196356) scale = 0.9943;
	if (run>=196357 && run<=198115) scale = 0.9930;
	if (run>=198116 && run<=198940) scale = 0.9956;
	if (run>=198941 && run<=199317) scale = 0.9931;
	if (run>=199318 && run<=199428) scale = 0.9930;
	if (run>=199429 && run<=199697) scale = 0.9936;
	if (run>=199698 && run<=199832) scale = 0.9920;
	if (run>=199833 && run<=199960) scale = 0.9946;
	if (run>=199961 && run<=200151) scale = 0.9927;
	if (run>=200152 && run<=200490) scale = 0.9947;
	if (run>=200491 && run<=200991) scale = 0.9927;
	if (run>=200992 && run<=201201) scale = 0.9920;
	if (run>=201202 && run<=201624) scale = 0.9918;
	if (run>=201625 && run<=201707) scale = 0.9917;
	if (run>=201708 && run<=202059) scale = 0.9937;
	if (run>=202060 && run<=202204) scale = 0.9923;
	if (run>=202205 && run<=202332) scale = 0.9946;
	if (run>=202333 && run<=202972) scale = 0.9922;
	if (run>=202973 && run<=203002) scale = 0.9934;
	if (run>=203003 && run<=203852) scale = 0.9781;
	if (run>=203853 && run<=204099) scale = 0.9928;
	if (run>=204100 && run<=204562) scale = 0.9928;
	if (run>=204563 && run<=205085) scale = 0.9940;
	if (run>=205086 && run<=205310) scale = 0.9924;
	if (run>=205311 && run<=205617) scale = 0.9925;
	if (run>=205618 && run<=205825) scale = 0.9928;
	if (run>=205826 && run<=206207) scale = 0.9929;
	if (run>=206208 && run<=206389) scale = 0.9909;
	if (run>=206390 && run<=206483) scale = 0.9950;
	if (run>=206484 && run<=206597) scale = 0.9934;
	if (run>=206598 && run<=206896) scale = 0.9923;
	if (run>=206897 && run<=207220) scale = 0.9939;
	if (run>=207221 && run<=207315) scale = 0.9934;
	if (run>=207316 && run<=207489) scale = 0.9930;
	if (run>=207490 && run<=207919) scale = 0.9933;
	if (run>=207920 && run<=208351) scale = 0.9928;
	if (run>=208352 && run<=208686) scale = 0.9928;
      } else if (fabs(scEta) >=2 && R9<0.94) {
	if (run>=190645 && run<=190781) scale = 1.0106;
	if (run>=190782 && run<=191042) scale = 1.0058;
	if (run>=191043 && run<=191720) scale = 1.0099;
	if (run>=191721 && run<=193833) scale = 1.0098;
	if (run>=193834 && run<=194116) scale = 1.0095;
	if (run>=194117 && run<=194427) scale = 1.0098;
	if (run>=194428 && run<=194618) scale = 1.0085;
	if (run>=194619 && run<=194789) scale = 1.0078;
	if (run>=194790 && run<=195111) scale = 1.0088;
	if (run>=195112 && run<=195377) scale = 1.0092;
	if (run>=195378 && run<=195398) scale = 1.0100;
	if (run>=195399 && run<=195657) scale = 1.0097;
	if (run>=195658 && run<=195918) scale = 1.0080;
	if (run>=195919 && run<=196198) scale = 1.0085;
	if (run>=196199 && run<=196356) scale = 1.0068;
	if (run>=196357 && run<=198115) scale = 1.0089;
	if (run>=198116 && run<=198940) scale = 1.0110;
	if (run>=198941 && run<=199317) scale = 1.0096;
	if (run>=199318 && run<=199428) scale = 1.0102;
	if (run>=199429 && run<=199697) scale = 1.0109;
	if (run>=199698 && run<=199832) scale = 1.0092;
	if (run>=199833 && run<=199960) scale = 1.0116;
	if (run>=199961 && run<=200151) scale = 1.0106;
	if (run>=200152 && run<=200490) scale = 1.0098;
	if (run>=200491 && run<=200991) scale = 1.0095;
	if (run>=200992 && run<=201201) scale = 1.0127;
	if (run>=201202 && run<=201624) scale = 1.0105;
	if (run>=201625 && run<=201707) scale = 1.0096;
	if (run>=201708 && run<=202059) scale = 1.0113;
	if (run>=202060 && run<=202204) scale = 1.0116;
	if (run>=202205 && run<=202332) scale = 1.0094;
	if (run>=202333 && run<=202972) scale = 1.0094;
	if (run>=202973 && run<=203002) scale = 1.0109;
	if (run>=203003 && run<=203852) scale = 1.0014;
	if (run>=203853 && run<=204099) scale = 1.0066;
	if (run>=204100 && run<=204562) scale = 1.0089;
	if (run>=204563 && run<=205085) scale = 1.0069;
	if (run>=205086 && run<=205310) scale = 1.0049;
	if (run>=205311 && run<=205617) scale = 1.0062;
	if (run>=205618 && run<=205825) scale = 1.0069;
	if (run>=205826 && run<=206207) scale = 1.0080;
	if (run>=206208 && run<=206389) scale = 1.0075;
	if (run>=206390 && run<=206483) scale = 1.0092;
	if (run>=206484 && run<=206597) scale = 1.0060;
	if (run>=206598 && run<=206896) scale = 1.0078;
	if (run>=206897 && run<=207220) scale = 1.0095;
	if (run>=207221 && run<=207315) scale = 1.0085;
	if (run>=207316 && run<=207489) scale = 1.0059;
	if (run>=207490 && run<=207919) scale = 1.0052;
	if (run>=207920 && run<=208351) scale = 1.0093;
	if (run>=208352 && run<=208686) scale = 1.0082;
      } else if (fabs(scEta) >=2 && R9>=0.94) {
	if (run>=190645 && run<=190781) scale = 0.9971;
	if (run>=190782 && run<=191042) scale = 0.9924;
	if (run>=191043 && run<=191720) scale = 0.9964;
	if (run>=191721 && run<=193833) scale = 0.9963;
	if (run>=193834 && run<=194116) scale = 0.9960;
	if (run>=194117 && run<=194427) scale = 0.9963;
	if (run>=194428 && run<=194618) scale = 0.9951;
	if (run>=194619 && run<=194789) scale = 0.9943;
	if (run>=194790 && run<=195111) scale = 0.9954;
	if (run>=195112 && run<=195377) scale = 0.9958;
	if (run>=195378 && run<=195398) scale = 0.9965;
	if (run>=195399 && run<=195657) scale = 0.9962;
	if (run>=195658 && run<=195918) scale = 0.9945;
	if (run>=195919 && run<=196198) scale = 0.9950;
	if (run>=196199 && run<=196356) scale = 0.9934;
	if (run>=196357 && run<=198115) scale = 0.9955;
	if (run>=198116 && run<=198940) scale = 0.9975;
	if (run>=198941 && run<=199317) scale = 0.9961;
	if (run>=199318 && run<=199428) scale = 0.9967;
	if (run>=199429 && run<=199697) scale = 0.9974;
	if (run>=199698 && run<=199832) scale = 0.9957;
	if (run>=199833 && run<=199960) scale = 0.9981;
	if (run>=199961 && run<=200151) scale = 0.9971;
	if (run>=200152 && run<=200490) scale = 0.9963;
	if (run>=200491 && run<=200991) scale = 0.9961;
	if (run>=200992 && run<=201201) scale = 0.9992;
	if (run>=201202 && run<=201624) scale = 0.9971;
	if (run>=201625 && run<=201707) scale = 0.9961;
	if (run>=201708 && run<=202059) scale = 0.9978;
	if (run>=202060 && run<=202204) scale = 0.9981;
	if (run>=202205 && run<=202332) scale = 0.9959;
	if (run>=202333 && run<=202972) scale = 0.9959;
	if (run>=202973 && run<=203002) scale = 0.9974;
	if (run>=203003 && run<=203852) scale = 0.9881;
	if (run>=203853 && run<=204099) scale = 0.9931;
	if (run>=204100 && run<=204562) scale = 0.9955;
	if (run>=204563 && run<=205085) scale = 0.9935;
	if (run>=205086 && run<=205310) scale = 0.9915;
	if (run>=205311 && run<=205617) scale = 0.9927;
	if (run>=205618 && run<=205825) scale = 0.9934;
	if (run>=205826 && run<=206207) scale = 0.9946;
	if (run>=206208 && run<=206389) scale = 0.9941;
	if (run>=206390 && run<=206483) scale = 0.9957;
	if (run>=206484 && run<=206597) scale = 0.9926;
	if (run>=206598 && run<=206896) scale = 0.9943;
	if (run>=206897 && run<=207220) scale = 0.9960;
	if (run>=207221 && run<=207315) scale = 0.9950;
	if (run>=207316 && run<=207489) scale = 0.9925;
	if (run>=207490 && run<=207919) scale = 0.9918;
	if (run>=207920 && run<=208351) scale = 0.9959;
	if (run>=208352 && run<=208686) scale = 0.9948;
      }
    } // end if dataset == Moriond2014 (jan22rereco)

    else if (targetDataset == "2011") {
      // values from http://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=176520
      double corr = 0;

      if ( fabs(scEta) < 1 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0014; 
	if (run>=170000 && run<=172619) corr = -0.0016;   
	if (run>=172620 && run<=173692) corr = -0.0017;  
	if (run>=175830 && run<=177139) corr = -0.0021;
	if (run>=177140 && run<=178421) corr = -0.0025;
	if (run>=178424 && run<=180252) corr = -0.0024;
      } else if (fabs(scEta) < 1 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0059; 
	if (run>=170000 && run<=172619) corr = 0.0046;   
	if (run>=172620 && run<=173692) corr = 0.0045;  
	if (run>=175830 && run<=177139) corr = 0.0042;
	if (run>=177140 && run<=178421) corr = 0.0038;
	if (run>=178424 && run<=180252) corr = 0.0039;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0045; 
	if (run>=170000 && run<=172619) corr = -0.0066;   
	if (run>=172620 && run<=173692) corr = -0.0058;  
	if (run>=175830 && run<=177139) corr = -0.0073;
	if (run>=177140 && run<=178421) corr = -0.0075;
	if (run>=178424 && run<=180252) corr = -0.0071;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0084; 
	if (run>=170000 && run<=172619) corr = 0.0063;   
	if (run>=172620 && run<=173692) corr = 0.0071;  
	if (run>=175830 && run<=177139) corr = 0.0056;
	if (run>=177140 && run<=178421) corr = 0.0054;
	if (run>=178424 && run<=180252) corr = 0.0058;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0082; 
	if (run>=170000 && run<=172619) corr = -0.0025;   
	if (run>=172620 && run<=173692) corr = -0.0035;  
	if (run>=175830 && run<=177139) corr = -0.0017;
	if (run>=177140 && run<=178421) corr = -0.0010;
	if (run>=178424 && run<=180252) corr = 0.0030;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = -0.0033; 
	if (run>=170000 && run<=172619) corr = 0.0024;   
	if (run>=172620 && run<=173692) corr = 0.0014;  
	if (run>=175830 && run<=177139) corr = 0.0032;
	if (run>=177140 && run<=178421) corr = 0.0040;
	if (run>=178424 && run<=180252) corr = 0.0079;
      } else if (fabs(scEta) >=2 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0064; 
	if (run>=170000 && run<=172619) corr = -0.0046;   
	if (run>=172620 && run<=173692) corr = -0.0029;  
	if (run>=175830 && run<=177139) corr = -0.0040;
	if (run>=177140 && run<=178421) corr = -0.0050;
	if (run>=178424 && run<=180252) corr = -0.0059;
      } else if (fabs(scEta) >= 2 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0042; 
	if (run>=170000 && run<=172619) corr = 0.0060;   
	if (run>=172620 && run<=173692) corr = 0.0077;  
	if (run>=175830 && run<=177139) corr = 0.0067;
	if (run>=177140 && run<=178421) corr = 0.0056;
	if (run>=178424 && run<=180252) corr = 0.0047;
      }
      
      scale = 1.0/(1.0+corr); //translate from corr to scale
    } //end if dataset == 2011 ( jan16 rereco )
    else {
      cout << "Warning: targetDataset " << targetDataset << " is not recogized. no corrections applied. \n";
    }

  } // end if type == 0

//Regression V0
  else if (type == 1) {
    if (targetDataset == "Moriond2013") {
      if ( fabs(scEta) < 1 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 1.0020;
        if (run>=190782 && run<=191042) scale = 1.0079;
        if (run>=191043 && run<=193555) scale = 0.9989;
        if (run>=193556 && run<=194150) scale = 0.9974;
        if (run>=194151 && run<=194532) scale = 0.9980;
        if (run>=194533 && run<=195113) scale = 0.9983;
        if (run>=195114 && run<=195915) scale = 0.9984;
        if (run>=195916 && run<=198115) scale = 0.9975;
        if (run>=198116 && run<=199803) scale = 1.0010;
        if (run>=199804 && run<=200048) scale = 1.0021;
        if (run>=200049 && run<=200151) scale = 1.0035;
        if (run>=200152 && run<=200490) scale = 1.0013;
        if (run>=200491 && run<=200531) scale = 1.0035;
        if (run>=200532 && run<=201656) scale = 1.0017;
        if (run>=201657 && run<=202305) scale = 1.0026;
        if (run>=202305 && run<=203002) scale = 1.0037;
        if (run>=203003 && run<=203984) scale = 1.0062;
        if (run>=203985 && run<=205085) scale = 1.0050;
        if (run>=205086 && run<=205310) scale = 1.0062;
        if (run>=205311 && run<=206207) scale = 1.0056;
        if (run>=206208 && run<=206483) scale = 1.0060;
        if (run>=206484 && run<=206597) scale = 1.0062;
        if (run>=206598 && run<=206896) scale = 1.0061;
        if (run>=206897 && run<=207220) scale = 1.0063;
        if (run>=207221 && run<=208686) scale = 1.0064;
      } else if (fabs(scEta) < 1 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9980;
        if (run>=190782 && run<=191042) scale = 1.0039;
        if (run>=191043 && run<=193555) scale = 0.9949;
        if (run>=193556 && run<=194150) scale = 0.9934;
        if (run>=194151 && run<=194532) scale = 0.9940;
        if (run>=194533 && run<=195113) scale = 0.9943;
        if (run>=195114 && run<=195915) scale = 0.9944;
        if (run>=195916 && run<=198115) scale = 0.9936;
        if (run>=198116 && run<=199803) scale = 0.9970;
        if (run>=199804 && run<=200048) scale = 0.9982;
        if (run>=200049 && run<=200151) scale = 0.9996;
        if (run>=200152 && run<=200490) scale = 0.9973;
        if (run>=200491 && run<=200531) scale = 0.9995;
        if (run>=200532 && run<=201656) scale = 0.9978;
        if (run>=201657 && run<=202305) scale = 0.9986;
        if (run>=202305 && run<=203002) scale = 0.9998;
        if (run>=203003 && run<=203984) scale = 1.0023;
        if (run>=203985 && run<=205085) scale = 1.0011;
        if (run>=205086 && run<=205310) scale = 1.0023;
        if (run>=205311 && run<=206207) scale = 1.0017;
        if (run>=206208 && run<=206483) scale = 1.0021;
        if (run>=206484 && run<=206597) scale = 1.0024;
        if (run>=206598 && run<=206896) scale = 1.0022;
        if (run>=206897 && run<=207220) scale = 1.0024;
        if (run>=207221 && run<=208686) scale = 1.0025;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 1.0032;
        if (run>=190782 && run<=191042) scale = 1.0063;
        if (run>=191043 && run<=193555) scale = 0.9998;
        if (run>=193556 && run<=194150) scale = 0.9954;
        if (run>=194151 && run<=194532) scale = 0.9965;
        if (run>=194533 && run<=195113) scale = 0.9984;
        if (run>=195114 && run<=195915) scale = 0.9977;
        if (run>=195916 && run<=198115) scale = 0.9965;
        if (run>=198116 && run<=199803) scale = 0.9999;
        if (run>=199804 && run<=200048) scale = 1.0008;
        if (run>=200049 && run<=200151) scale = 1.0017;
        if (run>=200152 && run<=200490) scale = 1.0003;
        if (run>=200491 && run<=200531) scale = 1.0017;
        if (run>=200532 && run<=201656) scale = 0.9999;
        if (run>=201657 && run<=202305) scale = 1.0003;
        if (run>=202305 && run<=203002) scale = 1.0010;
        if (run>=203003 && run<=203984) scale = 1.0022;
        if (run>=203985 && run<=205085) scale = 1.0045;
        if (run>=205086 && run<=205310) scale = 1.0046;
        if (run>=205311 && run<=206207) scale = 1.0034;
        if (run>=206208 && run<=206483) scale = 1.0036;
        if (run>=206484 && run<=206597) scale = 1.0033;
        if (run>=206598 && run<=206896) scale = 1.0021;
        if (run>=206897 && run<=207220) scale = 1.0033;
        if (run>=207221 && run<=208686) scale = 1.0036;     
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9919;
        if (run>=190782 && run<=191042) scale = 0.9951;
        if (run>=191043 && run<=193555) scale = 0.9885;
        if (run>=193556 && run<=194150) scale = 0.9841;
        if (run>=194151 && run<=194532) scale = 0.9852;
        if (run>=194533 && run<=195113) scale = 0.9872;
        if (run>=195114 && run<=195915) scale = 0.9864;
        if (run>=195916 && run<=198115) scale = 0.9852;
        if (run>=198116 && run<=199803) scale = 0.9886;
        if (run>=199804 && run<=200048) scale = 0.9895;
        if (run>=200049 && run<=200151) scale = 0.9905;
        if (run>=200152 && run<=200490) scale = 0.9890;
        if (run>=200491 && run<=200531) scale = 0.9905;
        if (run>=200532 && run<=201656) scale = 0.9887;
        if (run>=201657 && run<=202305) scale = 0.9891;
        if (run>=202305 && run<=203002) scale = 0.9897;
        if (run>=203003 && run<=203984) scale = 0.9914;
        if (run>=203985 && run<=205085) scale = 0.9938;
        if (run>=205086 && run<=205310) scale = 0.9939;
        if (run>=205311 && run<=206207) scale = 0.9926;
        if (run>=206208 && run<=206483) scale = 0.9929;
        if (run>=206484 && run<=206597) scale = 0.9925;
        if (run>=206598 && run<=206896) scale = 0.9913;
        if (run>=206897 && run<=207220) scale = 0.9926;
        if (run>=207221 && run<=208686) scale = 0.9929;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 0.9945;
        if (run>=190782 && run<=191042) scale = 0.9996;
        if (run>=191043 && run<=193555) scale = 0.9968;
        if (run>=193556 && run<=194150) scale = 0.9969;
        if (run>=194151 && run<=194532) scale = 0.9986;
        if (run>=194533 && run<=195113) scale = 1.0006;
        if (run>=195114 && run<=195915) scale = 1.0010;
        if (run>=195916 && run<=198115) scale = 1.0020;
        if (run>=198116 && run<=199803) scale = 0.9963;
        if (run>=199804 && run<=200048) scale = 0.9965;
        if (run>=200049 && run<=200151) scale = 0.9992;
        if (run>=200152 && run<=200490) scale = 0.9991;
        if (run>=200491 && run<=200531) scale = 0.9995;
        if (run>=200532 && run<=201656) scale = 0.9978;
        if (run>=201657 && run<=202305) scale = 0.9987;
        if (run>=202305 && run<=203002) scale = 1.0003;
        if (run>=203003 && run<=203984) scale = 0.9998;
        if (run>=203985 && run<=205085) scale = 1.0007;
        if (run>=205086 && run<=205310) scale = 1.0063;
        if (run>=205311 && run<=206207) scale = 1.0020;
        if (run>=206208 && run<=206483) scale = 0.9997;
        if (run>=206484 && run<=206597) scale = 1.0047;
        if (run>=206598 && run<=206896) scale = 0.9984;
        if (run>=206897 && run<=207220) scale = 0.9988;
        if (run>=207221 && run<=208686) scale = 1.0031;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9881;
        if (run>=190782 && run<=191042) scale = 0.9932;
        if (run>=191043 && run<=193555) scale = 0.9904;
        if (run>=193556 && run<=194150) scale = 0.9905;
        if (run>=194151 && run<=194532) scale = 0.9922;
        if (run>=194533 && run<=195113) scale = 0.9943;
        if (run>=195114 && run<=195915) scale = 0.9946;
        if (run>=195916 && run<=198115) scale = 0.9956;
        if (run>=198116 && run<=199803) scale = 0.9899;
        if (run>=199804 && run<=200048) scale = 0.9901;
        if (run>=200049 && run<=200151) scale = 0.9928;
        if (run>=200152 && run<=200490) scale = 0.9927;
        if (run>=200491 && run<=200531) scale = 0.9931;
        if (run>=200532 && run<=201656) scale = 0.9914;
        if (run>=201657 && run<=202305) scale = 0.9923;
        if (run>=202305 && run<=203002) scale = 0.9940;
        if (run>=203003 && run<=203984) scale = 0.9919;
        if (run>=203985 && run<=205085) scale = 0.9928;
        if (run>=205086 && run<=205310) scale = 0.9985;
        if (run>=205311 && run<=206207) scale = 0.9941;
        if (run>=206208 && run<=206483) scale = 0.9918;
        if (run>=206484 && run<=206597) scale = 0.9969;
        if (run>=206598 && run<=206896) scale = 0.9905;
        if (run>=206897 && run<=207220) scale = 0.9909;
        if (run>=207221 && run<=208686) scale = 0.9953;
      } else if (fabs(scEta) >=2 && R9<0.94) {
        if (run>=190645 && run<=190781) scale = 0.9965;
        if (run>=190782 && run<=191042) scale = 1.0010;
        if (run>=191043 && run<=193555) scale = 0.9987;
        if (run>=193556 && run<=194150) scale = 0.9988;
        if (run>=194151 && run<=194532) scale = 0.9994;
        if (run>=194533 && run<=195113) scale = 0.9999;
        if (run>=195114 && run<=195915) scale = 1.0004;
        if (run>=195916 && run<=198115) scale = 0.9992;
        if (run>=198116 && run<=199803) scale = 1.0044;
        if (run>=199804 && run<=200048) scale = 1.0060;
        if (run>=200049 && run<=200151) scale = 1.0101;
        if (run>=200152 && run<=200490) scale = 1.0073;
        if (run>=200491 && run<=200531) scale = 1.0106;
        if (run>=200532 && run<=201656) scale = 1.0069;
        if (run>=201657 && run<=202305) scale = 1.0121;
        if (run>=202305 && run<=203002) scale = 1.0144;
        if (run>=203003 && run<=203984) scale = 1.0128;
        if (run>=203985 && run<=205085) scale = 1.0114;
        if (run>=205086 && run<=205310) scale = 1.0213;
        if (run>=205311 && run<=206207) scale = 1.0162;
        if (run>=206208 && run<=206483) scale = 1.0186;
        if (run>=206484 && run<=206597) scale = 1.0177;
        if (run>=206598 && run<=206896) scale = 1.0173;
        if (run>=206897 && run<=207220) scale = 1.0185;
        if (run>=207221 && run<=208686) scale = 1.0195;
      } else if (fabs(scEta) >=2 && R9>=0.94) {
        if (run>=190645 && run<=190781) scale = 0.9862;
        if (run>=190782 && run<=191042) scale = 0.9907;
        if (run>=191043 && run<=193555) scale = 0.9884;
        if (run>=193556 && run<=194150) scale = 0.9885;
        if (run>=194151 && run<=194532) scale = 0.9891;
        if (run>=194533 && run<=195113) scale = 0.9896;
        if (run>=195114 && run<=195915) scale = 0.9900;
        if (run>=195916 && run<=198115) scale = 0.9889;
        if (run>=198116 && run<=199803) scale = 0.9941;
        if (run>=199804 && run<=200048) scale = 0.9957;
        if (run>=200049 && run<=200151) scale = 0.9999;
        if (run>=200152 && run<=200490) scale = 0.9970;
        if (run>=200491 && run<=200531) scale = 1.0004;
        if (run>=200532 && run<=201656) scale = 0.9967;
        if (run>=201657 && run<=202305) scale = 1.0018;
        if (run>=202305 && run<=203002) scale = 1.0042;
        if (run>=203003 && run<=203984) scale = 1.0027;
        if (run>=203985 && run<=205085) scale = 1.0013;
        if (run>=205086 && run<=205310) scale = 1.0113;
        if (run>=205311 && run<=206207) scale = 1.0061;
        if (run>=206208 && run<=206483) scale = 1.0086;
        if (run>=206484 && run<=206597) scale = 1.0077;
        if (run>=206598 && run<=206896) scale = 1.0073;
        if (run>=206897 && run<=207220) scale = 1.0085;
        if (run>=207221 && run<=208686) scale = 1.0094;
      }
    } // end if dataset == Moriond2013 (j13,a06,a24 rereco, r12c pr-v2, r12d pr-v1)

    else if (targetDataset == "2011") {

      double corr = 0;
      
      if ( fabs(scEta) < 1 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0014; 
	if (run>=170000 && run<=172619) corr = -0.0016;   
	if (run>=172620 && run<=173692) corr = -0.0017;  
	if (run>=175830 && run<=177139) corr = -0.0021;
	if (run>=177140 && run<=178421) corr = -0.0025;
	if (run>=178424 && run<=180252) corr = -0.0024;
      } else if (fabs(scEta) < 1 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0059; 
	if (run>=170000 && run<=172619) corr = 0.0046;   
	if (run>=172620 && run<=173692) corr = 0.0045;  
	if (run>=175830 && run<=177139) corr = 0.0042;
	if (run>=177140 && run<=178421) corr = 0.0038;
	if (run>=178424 && run<=180252) corr = 0.0039;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0045; 
	if (run>=170000 && run<=172619) corr = -0.0066;   
	if (run>=172620 && run<=173692) corr = -0.0058;  
	if (run>=175830 && run<=177139) corr = -0.0073;
	if (run>=177140 && run<=178421) corr = -0.0075;
	if (run>=178424 && run<=180252) corr = -0.0071;
      } else if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0084; 
	if (run>=170000 && run<=172619) corr = 0.0063;   
	if (run>=172620 && run<=173692) corr = 0.0071;  
	if (run>=175830 && run<=177139) corr = 0.0056;
	if (run>=177140 && run<=178421) corr = 0.0054;
	if (run>=178424 && run<=180252) corr = 0.0058;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2  && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0082; 
	if (run>=170000 && run<=172619) corr = -0.0025;   
	if (run>=172620 && run<=173692) corr = -0.0035;  
	if (run>=175830 && run<=177139) corr = -0.0017;
	if (run>=177140 && run<=178421) corr = -0.0010;
	if (run>=178424 && run<=180252) corr = 0.0030;
      } else if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = -0.0033; 
	if (run>=170000 && run<=172619) corr = 0.0024;   
	if (run>=172620 && run<=173692) corr = 0.0014;  
	if (run>=175830 && run<=177139) corr = 0.0032;
	if (run>=177140 && run<=178421) corr = 0.0040;
	if (run>=178424 && run<=180252) corr = 0.0079;
      } else if (fabs(scEta) >=2 && R9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0064; 
	if (run>=170000 && run<=172619) corr = -0.0046;   
	if (run>=172620 && run<=173692) corr = -0.0029;  
	if (run>=175830 && run<=177139) corr = -0.0040;
	if (run>=177140 && run<=178421) corr = -0.0050;
	if (run>=178424 && run<=180252) corr = -0.0059;
      } else if (fabs(scEta) >= 2 && R9>=0.94) {
        if (run>=160431 && run<=167913) corr = 0.0042; 
	if (run>=170000 && run<=172619) corr = 0.0060;   
	if (run>=172620 && run<=173692) corr = 0.0077;  
	if (run>=175830 && run<=177139) corr = 0.0067;
	if (run>=177140 && run<=178421) corr = 0.0056;
	if (run>=178424 && run<=180252) corr = 0.0047;
      }
      
      scale = 1.0/(1.0+corr); //translate from corr to scale
    }
    else {
      cout << "Warning: targetDataset " << targetDataset << " is not recogized. no corrections applied. \n";
    }

  }

  double dsigMC = 0.0;

  if (type == 0) {
    if (targetDataset == "HCP2012") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0099;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0103;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0219;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0158;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0222;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0298;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0318;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0302;        
    }
    else if (targetDataset == "ICHEP2012") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0119;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0107;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0240;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0149;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0330;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0375;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0602;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0607;        
    }
    else if (targetDataset == "Moriond2014") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0114;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0092;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0246;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0133;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0229;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0200;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0236;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0246;        
    }
    else if (targetDataset == "2011") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0096;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0074;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0196;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0141;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0279;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0268;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0301;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0293;        
    }
    else {
      cout << "Warning: targetDataset " << targetDataset << " is not recogized. no corrections applied. \n";
    }
  }
  else if (type == 1) {
    if (targetDataset == "HCP2012") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0103;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0090;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0190;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0156;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0269;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0287;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0364;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0321;        
    } 
    else if (targetDataset == "Moriond2013") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0112;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0100;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0208;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0162;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0288;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0285;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0368;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0360;        
    }
    else if (targetDataset == "2011") {
      if (fabs(scEta) <1 && R9<0.94) dsigMC = 0.0096;
      if (fabs(scEta) <1 && R9>=0.94) dsigMC = 0.0074;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9<0.94) dsigMC = 0.0196;
      if (fabs(scEta) < 1.479 && fabs(scEta) >= 1 && R9>=0.94) dsigMC = 0.0141;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9<0.94) dsigMC = 0.0279;
      if (fabs(scEta) >= 1.479 && fabs(scEta) < 2 && R9>=0.94) dsigMC = 0.0268;
      if (fabs(scEta) >=2 && R9<0.94) dsigMC = 0.0301;
      if (fabs(scEta) >=2 && R9>=0.94) dsigMC = 0.0293;        
    }
    else {
      cout << "Warning: targetDataset " << targetDataset << " is not recogized. no corrections applied. \n";
    }
  }

  //perform corrections
  double newEnergy = inputEnergy;

  //do scale corr
  if (!isMC) {
    newEnergy = inputEnergy * scale;
  }

  //smearing
  if (isMC) {
    newEnergy = newEnergy * random->Gaus(1,dsigMC);
  }
    
  return newEnergy;

}

#endif
