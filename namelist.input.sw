&sw_model
   config_test_case = 5
   config_time_integration = 'RK4'
   config_dt = 172.8
   config_ntimesteps = 7500
   config_output_interval = 500
   config_stats_interval = 0
   config_h_mom_eddy_visc2  = 0.0
   config_h_mom_eddy_visc4  = 0.0
   config_h_tracer_eddy_diff2  = 0.0
   config_h_tracer_eddy_diff4  = 0.0
   config_thickness_adv_order = 2
   config_tracer_adv_order = 2
   config_positive_definite = .false.
   config_monotonic = .false.
/

&io
   config_input_name = 'grid.nc'
   config_output_name = 'output.nc'
   config_restart_name = 'restart.nc'
/

&restart
   config_restart_interval = 3000
   config_do_restart = .false.
   config_restart_time = 1036800.0
/
