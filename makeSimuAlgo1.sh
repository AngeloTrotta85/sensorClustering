#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1

rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"
rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"

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
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			$EXEC -f inputScenario/inputTestAlgo1.dat -g inputScenario/inputTestAlgo1.dat -s 100 -e 0 -n $sensors -l 1 &>/dev/null
		
			NOW_T=`date +"%F %T"`
			echo "$NOW_T s:$sensors l:$lambda r:$runs"
		
			ALGO1=`$EXEC -f inputScenario/inputTestAlgo1.dat -e 4 -l $lambda | grep StatMaxCorr`
			ALGO1CORR=`echo $ALGO1 | awk '{printf $2}'`
			ALGO1DIST=`echo $ALGO1 | awk '{printf $3}'`
			echo "$ALGO1CORR" >> "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
			echo "$ALGO1DIST" >> "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
			#echo ""
		done
		
		ALGO1CORR_OK=`cat "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		ALGO1DIST_OK=`cat "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"					
		
		echo "$lambda $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_s${sensors}.data"
		echo "$sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_l${lambda}.data"
		echo "$lambda $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_s${sensors}.data"
		echo "$sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_l${lambda}.data"
		
		echo "$lambda $sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
		echo "$lambda $sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
	done
	
	echo "" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

