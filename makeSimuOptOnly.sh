#! /bin/bash

OUTPUT_DIR="statsOpt"
EXEC="Release/SensorsClustering"
N_RUNS=$1
MAX_OPT_OP=$2

rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"


for lambda in $3
#for lambda in {$3..25..$4}
#for (( lambda=MIN_SENSORS; lambda<=12; lambda+=lambda ))
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
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			SEED=`echo $(( runs + (lambda * 1000) + (sensors * 100000) ))`
			#echo "$NOW_T s:$sensors l:$lambda r:$runs  ->  SEED = $SEED"
		
			NOW_T=`date +"%F %T"`
			echo -n "$NOW_T s:$sensors l:$lambda r:$runs SEED:$SEED  ->  "
			
			
			$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -g inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -s 200 -e 0 -n $sensors -z $SEED -l 1 &>/dev/null
		
			echo -n "Algo1... "
						
			echo -n "OK - Opt... "
			OPT=`$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -t $MAX_OPT_OP -e 5 -l $lambda | grep StatMaxCorr`
			OPTCORR=`echo $OPT | awk '{printf $2}'`
			OPTDIST=`echo $OPT | awk '{printf $3}'`
			echo "$OPTCORR" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			echo "$OPTDIST" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			
			#EXPECTED_OP=`echo $(( (sensors/lambda) ** sensors ))`
			#if [ $EXPECTED_OP -lt $MAX_OPT_OP ]
			#then
			#	echo -n "OK - Opt ($EXPECTED_OP expected operation)... "
			#	OPT=`$EXEC -f inputScenario/inputTest.dat -e 5 -l $lambda | grep StatMaxCorr`
			#	OPTCORR=`echo $OPT | awk '{printf $2}'`
			#	OPTDIST=`echo $OPT | awk '{printf $3}'`
			#	echo "$OPTCORR" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			#	echo "$OPTDIST" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			#else
			#	echo "0" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			#	echo "0" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			#	echo -n "OK - Skipping Opt ($EXPECTED_OP expected operation) "
			#fi
		
			
		
			echo -n "OK - Random... "

			
		
			echo -n "OK - Salsiccia... "

		
			echo "OK"
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

