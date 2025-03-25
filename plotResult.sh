#! /bin/sh

########################
STA=B17L
PS=Joint-${STA}.ps
Depth1=0
Depth2=100
################### jiahui #######
gmt gmtset MAP_FRAME_TYPE fancy
gmt gmtset MAP_FRAME_PEN 1p
gmt gmtset MAP_FRAME_WIDTH 2p
gmt gmtset FONT_TITLE 10p
gmt gmtset FONT_ANNOT 7.0
gmt gmtset FONT_LABEL 8.5
gmt gmtset MAP_TITLE_OFFSET 0.2c
gmt gmtset MAP_LABEL_OFFSET 0.1c
################## jiahui #######

gmt makecpt -Cseis -T0/0.5/0.05 -Z -I > pdf.cpt

gmt psbasemap -JX1.5i/1.5i -R5/55/0.6/4.5 -Ba5:"Period (sec)":/a0.5f0.1:"Phase and Group velocity (km/s)":nSWe -K -P > ${PS}
cat ./SWD_phase.obs| awk '{print $1,$2}' | gmt psxy -J -R -Sc0.065i -W0.5p -K -O >> ${PS}
cat ./OUT/data_best_phase.out | awk '{print $1,$2}' | gmt psxy -J -R -W1p,red -K -O >> ${PS}
cat ./SWD_group.obs| awk '{print $1,$2}' | gmt psxy -J -R -Sc0.065i -W0.5p -K -O >> ${PS}
cat ./OUT/data_best_group.out | awk '{print $1,$2}' | gmt psxy -J -R -W1p,red -K -O >> ${PS}

gmt psxy -T -JX1.5i/1i -R-5/25/-0.8/1.2 -Ba5f1:"Time after P (s)":/a0S -K -O -Y2.0i  >>  ${PS}
cat ./RF2.obs | awk '{if(NR>1) print $1,$2*2}' | gmt psxy -J -R -W1p -K -O >> ${PS}
cat ./OUT/data_bestrg2.out | awk '{print $1,$2*2}' | gmt psxy -J -R -W1p,red -K -O  >>  ${PS}
cat ./RF3.obs | awk '{if(NR>1) print $1,$2*2}' | gmt psxy -J -R -W1p -K -O >> ${PS}
cat ./OUT/data_bestrg3.out | awk '{print $1,$2*2}' | gmt psxy -J -R -W1p,red -K -O  >>  ${PS}


gmt psbasemap -JX1.5i/-2.5i -R2.0/5.0/${Depth1}/${Depth2} -Bxa0.5f0.25+l"Vs (km/s)" -Bya10f2+l"Depth (km)" -BnSwE -K -O -Y-2.15i -X2i>>  ${PS} #8.2i

IN=./OUT/posteriorg.out
max=`cat ${IN} | awk '{if(NR>1 && $2<=100) print $3}' | sort -g | tail -1`
min=`cat ${IN} | awk '{if(NR>1 && $2<=100) print $3}' | sort -g | head -1`
cat ${IN} | awk '{if(NR>1) print $1,$2,($3-'${min}')/('${max}'-'${min}')}' > PDF.norm
gmt surface -Ll0.0 -Lu0.5 -T0.75 -R -I0.02/0.1 PDF.norm -GPDF.grd
gmt grdimage PDF.grd  -J -R -Cpdf.cpt -K -O >> ${PS}

cat ./OUT/Averageg.out | awk '{print $2,$1}' | gmt psxy -J -R -W2p,red -K -O >>  ${PS}
cat ./inputmodel.txt | awk '{print $2,$1}' | gmt psxy -J -R -W2p,black -K -O >>  ${PS}

gmt psscale -J -R -Cpdf.cpt -DjTR+w0c/0c+o0.5c/-0.45c+m+h -Ba0.1+l"PDF" -O >> ${PS} #-K

gmt psconvert ${PS} -Tf
rm -rf ${PS}
rm -rf PDF*
rm -rf pdf.cpt