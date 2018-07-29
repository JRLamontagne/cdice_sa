#!/bin/bash

cat CDICE_output_*txt>a.txt
cut -d' ' -f3 a.txt > util1.txt
cut -d' ' -f4 a.txt > pop1.txt
cut -d' ' -f5 a.txt > util_disc1.txt
cut -d' ' -f6 a.txt > real_intr1.txt
cut -d' ' -f7 a.txt > NPV_Dam1.txt
cut -d' ' -f8 a.txt > NPV_Abat1.txt
cut -d' ' -f9 a.txt > Abat1.txt
cut -d' ' -f10 a.txt > Dam1.txt
cut -d' ' -f11 a.txt > GWP1.txt
cut -d' ' -f12 a.txt > SCC1.txt
cut -d' ' -f13 a.txt > Tatm1.txt
cut -d' ' -f14 a.txt > Mat1.txt
cut -d' ' -f15 a.txt > forc1.txt
cut -d' ' -f16 a.txt > caa1.txt
cut -d' ' -f17 a.txt > caa_up_step1.txt
cut -d' ' -f18 a.txt > ygross1.txt

awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' util1.txt > util.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' pop1.txt > pop.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' util_disc1.txt > util_disc.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' real_intr1.txt > real_intr.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' NPV_Dam1.txt > NPV_Dam.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' NPV_Abat1.txt > NPV_Abat.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' Abat1.txt > Abat.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' Dam1.txt > Dam.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' GWP1.txt > GWP.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' SCC1.txt > SCC.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' Tatm1.txt > Tatm.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' Mat1.txt > Mat.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' forc1.txt > forc.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' caa1.txt > caa.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' caa_up_step1.txt > caa_up_step.txt
awk '{if(n==99){n=0;print $0}else{printf "%s ",$0;n++}}' ygross1.txt > ygross.txt

rm util1.txt
rm pop1.txt
rm util_disc1.txt
rm real_intr1.txt
rm NPV_Dam1.txt
rm NPV_Abat1.txt
rm Abat1.txt
rm Dam1.txt
rm GWP1.txt
rm SCC1.txt
rm Tatm1.txt
rm Mat1.txt
rm forc1.txt
rm caa1.txt
rm caa_up_step1.txt
rm ygross1.txt
rm a.txt