#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1
MAX_OPT_OP=$2

rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"
rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"

rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"

rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"

rm -rf "$OUTPUT_DIR/salscorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/salscorr_OK_3D.data"


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
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/salscorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/salsdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			SEED=`echo $(( runs + (lambda * 1000) + (sensors * 100000) ))`
			#echo "$NOW_T s:$sensors l:$lambda r:$runs  ->  SEED = $SEED"
		
			NOW_T=`date +"%F %T"`
			echo -n "$NOW_T s:$sensors l:$lambda r:$runs SEED:$SEED  ->  "
			
			
			$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -g inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -s 200 -e 0 -n $sensors -z $SEED -l 1 &>/dev/null
		
			echo -n "Algo1... "
			ALGO1=`$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -e 4 -l $lambda | grep StatMaxCorr`
			ALGO1CORR=`echo $ALGO1 | awk '{printf $2}'`
			ALGO1DIST=`echo $ALGO1 | awk '{printf $3}'`
			echo "$ALGO1CORR" >> "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
			echo "$ALGO1DIST" >> "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
			#echo "Algo1 corr:$ALGO1CORR dist:$ALGO1DIST"
						
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
			RAND=`$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -e 7 -l $lambda | grep StatMaxCorr`
			RANDCORR=`echo $RAND | awk '{printf $2}'`
			RANDDIST=`echo $RAND | awk '{printf $3}'`
			echo "$RANDCORR" >> "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
			echo "$RANDDIST" >> "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
			
		
			echo -n "OK - Salsiccia... "
			SALS=`$EXEC -f inputScenario/inputTest_l${lambda}_s${sensors}_r${runs}.dat -e 2 -l $lambda -i 100 | grep StatMaxCorr`
			SALSCORR=`echo $SALS | awk '{printf $2}'`
			SALSDIST=`echo $SALS | awk '{printf $3}'`
			echo "$SALSCORR" >> "$OUTPUT_DIR/salscorr_${lambda}_${sensors}.data"
			echo "$SALSDIST" >> "$OUTPUT_DIR/salsdist_${lambda}_${sensors}.data"
		
			echo "OK"
		done
		
		ALGO1CORR_OK=`cat "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		ALGO1DIST_OK=`cat "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		OPTCORR_OK=`cat "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		OPTDIST_OK=`cat "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		RANDCORR_OK=`cat "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		RANDDIST_OK=`cat "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		SALSCORR_OK=`cat "$OUTPUT_DIR/salscorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		SALSDIST_OK=`cat "$OUTPUT_DIR/salsdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/salscorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/salsdist_${lambda}_${sensors}.data"
				
		
		
		echo "$lambda $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_s${sensors}.data"
		echo "$sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_l${lambda}.data"
		echo "$lambda $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_s${sensors}.data"
		echo "$sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_l${lambda}.data"
		
		echo "$lambda $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_s${sensors}.data"
		echo "$sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_l${lambda}.data"
		echo "$lambda $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_s${sensors}.data"
		echo "$sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_l${lambda}.data"

		echo "$lambda $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_s${sensors}.data"
		echo "$sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_l${lambda}.data"
		echo "$lambda $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_s${sensors}.data"
		echo "$sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_l${lambda}.data"

		echo "$lambda $SALSCORR_OK" >> "$OUTPUT_DIR/salscorr_OK_s${sensors}.data"
		echo "$sensors $SALSCORR_OK" >> "$OUTPUT_DIR/salscorr_OK_l${lambda}.data"
		echo "$lambda $SALSDIST_OK" >> "$OUTPUT_DIR/salsdist_OK_s${sensors}.data"
		echo "$sensors $SALSDIST_OK" >> "$OUTPUT_DIR/salsdist_OK_l${lambda}.data"
		
		echo "$lambda $sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
		echo "$lambda $sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
		
		echo "$lambda $sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
		echo "$lambda $sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_3D.data"
		
		echo "$lambda $sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
		echo "$lambda $sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_3D.data"
		
		echo "$lambda $sensors $SALSCORR_OK" >> "$OUTPUT_DIR/salscorr_OK_3D.data"
		echo "$lambda $sensors $SALSDIST_OK" >> "$OUTPUT_DIR/salsdist_OK_3D.data"
	done
	
	echo "" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
	
	echo "" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/randdist_OK_3D.data"
	
	echo "" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/optdist_OK_3D.data"
	
	echo "" >> "$OUTPUT_DIR/salscorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/salsdist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

