library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity snn_architecture is
    port(clk        : in  std_logic;                        -- input clock
         reset_n     : in  std_logic;                        -- reset
         -- spike train input
         sp_0        : in  std_logic;                        -- spike signal for feature 0
         sp_1        : in  std_logic;                        -- spike signal for feature 1
         sp_2        : in  std_logic;                        -- spike signal for feature 2
         sp_3        : in  std_logic;                        -- spike signal for feature 3
         sp_4        : in  std_logic;                        -- spike signal for feature 4
         sp_5        : in  std_logic;                        -- spike signal for feature 5
         sp_6        : in  std_logic;                        -- spike signal for feature 6
         sp_7        : in  std_logic;                        -- spike signal for feature 7
         sp_8        : in  std_logic;                        -- spike signal for feature 8
         sp_9        : in  std_logic;                        -- spike signal for feature 9
         -- output signal
         spike_out_0   : out std_logic;                         -- output spike signal for neuron 0
         spike_out_1   : out std_logic;                          -- output spike signal for neuron 1
         voltage_out_0 : out integer;
         voltage_out_1 : out integer
    );
end snn_architecture;

architecture behave of snn_architecture is
    -- internal signals for neurons
    signal output_spike_0 : std_logic;                       -- output spike signal for neuron 0
    signal output_spike_1 : std_logic;                       -- output spike signal for neuron 1
    signal v_out_0 : integer   :=   0;
    signal v_out_1 : integer   :=   0;
    constant v_th : integer := 100; -- voltage threshold for spike generation
begin

    -- Output neuron 0
    output_neuron_0: entity work.neuron
        generic map(w_0   => 10,
                    w_1   => 8,
                    w_2   => 5,
                    w_3   => 20,
                    w_4   => 25,
                    w_5   => 10,
                    w_6   => 15,
                    w_7   => 8,
                    w_8   => 7,
                    w_9   => 5,
                    bias  => 0,
                    v_th  => v_th)
        port map(clk        => clk,
                 reset      => reset_n,
                 sp_0       => sp_0,
                 sp_1       => sp_1,
                 sp_2       => sp_2,
                 sp_3       => sp_3,
                 sp_4       => sp_4,
                 sp_5       => sp_5,
                 sp_6       => sp_6,
                 sp_7       => sp_7,
                 sp_8       => sp_8,
                 sp_9       => sp_9,
                 neuron_reset => '0',
                 spike_out  => output_spike_0,
                 voltage_out => v_out_0);

    -- Output neuron 1
    output_neuron_1: entity work.neuron
        generic map(w_0   => 20,
                    w_1   => 15,
                    w_2   => 5,
                    w_3   => 7,
                    w_4   => 12,
                    w_5   => 13,
                    w_6   => 8,
                    w_7   => 5,
                    w_8   => 15,
                    w_9   => 6,
                    bias  => 0,
                    v_th  => v_th)
        port map(clk        => clk,
                 reset      => reset_n,
                 sp_0       => sp_0,
                 sp_1       => sp_1,
                 sp_2       => sp_2,
                 sp_3       => sp_3,
                 sp_4       => sp_4,
                 sp_5       => sp_5,
                 sp_6       => sp_6,
                 sp_7       => sp_7,
                 sp_8       => sp_8,
                 sp_9       => sp_9,
                 neuron_reset => '0',
                 spike_out  => output_spike_1,
                 voltage_out => v_out_1);

    -- Output spike signals
    spike_out_0   <= output_spike_0;
    spike_out_1   <= output_spike_1;
    voltage_out_0 <= v_out_0;
    voltage_out_1 <= v_out_1;
    
end behave;