#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1
MAX_OPT_OP=$2


rm -rf "$OUTPUT_DIR/rrcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/rrcorr_OK_3D.data"


#for lambda in 3 4 5
for lambda in {3..25}
do
	MIN_SENSORS=`echo $(( lambda * 2 ))`
	#echo "MIN_SENSORS $MIN_SENSORS"
	
	#for sensors in {0..150..2}
	#for (( sensors=MIN_SENSORS; sensors<=120; sensors+=2 ))
	for (( sensors=MIN_SENSORS; sensors<=150; sensors+=lambda ))
	do
		if [ $sensors -lt $MIN_SENSORS ]
		then
			continue
		fi
		#echo "s:$sensors l:$lambda"
		#continue
		
		
		rm -rf "$OUTPUT_DIR/rrcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/rrdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			SEED=`echo $(( runs + (lambda * 1000) + (sensors * 100000) ))`
			#echo "$NOW_T s:$sensors l:$lambda r:$runs  ->  SEED = $SEED"
		
			NOW_T=`date +"%F %T"`
			echo -n "$NOW_T s:$sensors l:$lambda r:$runs SEED:$SEED  ->  "
			
			SCENARIO_FN="inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat"
			
			if [ ! -f "$SCENARIO_FN" ]
			then
				$EXEC -f ${SCENARIO_FN} -g ${SCENARIO_FN} -s 200 -e 0 -n $sensors -z $SEED -l 1 &>/dev/null
			fi
			
		
			echo -n "Algo1... "
						
			echo -n "OK - Opt... "
			
		
			echo -n "OK - Random... "
			
		
			echo -n "OK - Salsiccia... "
			
		
			echo -n "OK - Kmean... "
			
			
			echo -n "OK - roundRobin... "
			RR=`$EXEC -f ${SCENARIO_FN} -e 9 -l $lambda | grep StatMaxCorr`
			RRCORR=`echo $RR | awk '{printf $2}'`
			RRDIST=`echo $RR | awk '{printf $3}'`
			echo "$RRCORR" >> "$OUTPUT_DIR/rrcorr_${lambda}_${sensors}.data"
			echo "$RRDIST" >> "$OUTPUT_DIR/rrdist_${lambda}_${sensors}.data"
			
		
			echo "OK"
		done
		
		
		RRCORR_OK=`cat "$OUTPUT_DIR/rrcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		RRDIST_OK=`cat "$OUTPUT_DIR/rrdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		
		rm -rf "$OUTPUT_DIR/rrcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/rrdist_${lambda}_${sensors}.data"

		echo "$lambda $RRCORR_OK" >> "$OUTPUT_DIR/rrcorr_OK_s${sensors}.data"
		echo "$sensors $RRCORR_OK" >> "$OUTPUT_DIR/rrcorr_OK_l${lambda}.data"
		echo "$lambda $RRDIST_OK" >> "$OUTPUT_DIR/rrdist_OK_s${sensors}.data"
		echo "$sensors $RRDIST_OK" >> "$OUTPUT_DIR/rrdist_OK_l${lambda}.data"
		
		
		
		echo "$lambda $sensors $RRCORR_OK" >> "$OUTPUT_DIR/rrcorr_OK_3D.data"
		echo "$lambda $sensors $RRDIST_OK" >> "$OUTPUT_DIR/rrdist_OK_3D.data"
	done
	
	
	echo "" >> "$OUTPUT_DIR/rrcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/rrdist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

