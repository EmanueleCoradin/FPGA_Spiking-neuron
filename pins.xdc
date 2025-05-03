# ----------------------- CLOCK -----------------------
set_property -dict {PACKAGE_PIN E3 IOSTANDARD LVCMOS33} [get_ports { CLK100MHZ }];
create_clock -add -period 10 -name sys_clk_pin -waveform {0 5} [get_ports { CLK100MHZ }];  # Create clock with 100MHz (10ns period)
#set_property CLOCK_DEDICATED_ROUTE FALSE [get_nets clock_IBUF];

# ----------------------- UART -----------------------
set_property -dict {PACKAGE_PIN A9  IOSTANDARD LVCMOS33} [get_ports { uart_txd_in }];

# set_property PACKAGE_PIN D10 [get_ports uart_rxd_out]
# set_property IOSTANDARD LVCMOS33 [get_ports uart_rxd_out]

# ----------------------- LEDs -----------------------
set_property -dict {PACKAGE_PIN H5     IOSTANDARD LVCMOS33 }  [get_ports {led[0]}];
set_property -dict {PACKAGE_PIN J5     IOSTANDARD LVCMOS33 }  [get_ports {led[1]}];
set_property -dict {PACKAGE_PIN T9     IOSTANDARD LVCMOS33 }  [get_ports {led[2]}];
set_property -dict {PACKAGE_PIN T10    IOSTANDARD LVCMOS33 }  [get_ports {led[3]}];
set_property -dict { PACKAGE_PIN G6    IOSTANDARD LVCMOS33 } [get_ports  { led[4] }]; #IO_L19P_T3_35 Sch=led0_r
set_property -dict { PACKAGE_PIN G3    IOSTANDARD LVCMOS33 } [get_ports  { led[5] }]; #IO_L20N_T3_35 Sch=led1_r
set_property -dict { PACKAGE_PIN J3    IOSTANDARD LVCMOS33 } [get_ports  { led[6] }]; #IO_L22P_T3_35 Sch=led2_r
set_property -dict { PACKAGE_PIN K1    IOSTANDARD LVCMOS33 } [get_ports  { led[7] }]; #IO_L23N_T3_35 Sch=led3_r
