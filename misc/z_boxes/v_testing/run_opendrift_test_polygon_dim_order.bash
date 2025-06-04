
loggerLevel="DEBUG"
jobRunNum=0
jobRunString="BLAH"

configFile="/data03/blaughli/tracking_output/Mercator_hindcast_1993_2018_2D_dt_120_kick_0p3___mercator_reanalysis12_1993_2018/z_config_files/configFile_000.config.yaml"

python /home/blaughli/tracking_project_v2/run_scripts/opendrift_run.py --configfile $configFile --jobrunnumber $jobRunNum --level $loggerLevel --jobrunstring $jobRunString 
