#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1

rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"

#for lambda in 3 4 5
for lambda in {3..20}
do
	MIN_SENSORS=`echo $(( lambda * 2 ))`
	#echo "MIN_SENSORS $MIN_SENSORS"
	
	for sensors in {0..120..2}
	#for (( sensors=MIN_SENSORS; sensors<=120; sensors+=2 ))
	do
		if [ $sensors -lt $MIN_SENSORS ]
		then
			continue
		fi
		#echo "s:$sensors l:$lambda"
		#continue
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			$EXEC -f inputScenario/inputTestRand.dat -g inputScenario/inputTestRand.dat -s 100 -e 0 -n $sensors -l 1 &>/dev/null
		
			NOW_T=`date +"%F %T"`
			echo "$NOW_T s:$sensors l:$lambda r:$runs"
		
			RAND=`$EXEC -f inputScenario/inputTestRand.dat -e 7 -l $lambda | grep StatMaxCorr`
			RANDCORR=`echo $RAND | awk '{printf $2}'`
			RANDDIST=`echo $RAND | awk '{printf $3}'`
			echo "$RANDCORR" >> "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
			echo "$RANDDIST" >> "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
			#echo ""
		done
		
		RANDCORR_OK=`cat "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		RANDDIST_OK=`cat "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"

		echo "$lambda $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_s${sensors}.data"
		echo "$sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_l${lambda}.data"
		echo "$lambda $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_s${sensors}.data"
		echo "$sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_l${lambda}.data"
		
		echo "$lambda $sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
		echo "$lambda $sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_3D.data"
	done
	
	echo "" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/randdist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

