#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1

rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"

#for lambda in 3 4 5
for lambda in {3..5}
do
	MIN_SENSORS=`echo $(( lambda * 2 ))`
	#echo "MIN_SENSORS $MIN_SENSORS"
	
	for sensors in {0..12..2}
	#for (( sensors=MIN_SENSORS; sensors<=120; sensors+=2 ))
	do
		if [ $sensors -lt $MIN_SENSORS ]
		then
			continue
		fi
		#echo "s:$sensors l:$lambda"
		#continue
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			$EXEC -f inputScenario/inputTestOpt.dat -g inputScenario/inputTestOpt.dat -s 100 -e 0 -n $sensors -l 1 &>/dev/null
		
			NOW_T=`date +"%F %T"`
			echo "$NOW_T s:$sensors l:$lambda r:$runs"
		
			OPT=`$EXEC -f inputScenario/inputTestOpt.dat -e 5 -l $lambda | grep StatMaxCorr`
			OPTCORR=`echo $OPT | awk '{printf $2}'`
			OPTDIST=`echo $OPT | awk '{printf $3}'`
			echo "$OPTCORR" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			echo "$OPTDIST" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
			#echo ""
		done
		
		OPTCORR_OK=`cat "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		OPTDIST_OK=`cat "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
				
		echo "$lambda $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_s${sensors}.data"
		echo "$sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_l${lambda}.data"
		echo "$lambda $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_s${sensors}.data"
		echo "$sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_l${lambda}.data"
		
		echo "$lambda $sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
		echo "$lambda $sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_3D.data"
	done
	
	echo "" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/optdist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

