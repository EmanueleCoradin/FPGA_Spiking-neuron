library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.STD_LOGIC_ARITH.ALL;
use IEEE.STD_LOGIC_UNSIGNED.ALL;

entity tb_snn_spikes is
end tb_snn_spikes;

architecture behavior of tb_snn_spikes is

    -- Component Declaration for the Unit Under Test (UUT)
    component snn_spikes
        port(
            clk         : in  std_logic;
            reset_n     : in  std_logic;
            sp_0        : in  std_logic;
            sp_1        : in  std_logic;
            sp_2        : in  std_logic;
            sp_3        : in  std_logic;
            sp_4        : in  std_logic;
            sp_5        : in  std_logic;
            sp_6        : in  std_logic;
            sp_7        : in  std_logic;
            sp_8        : in  std_logic;
            sp_9        : in  std_logic;
            spike_out_0 : out std_logic;
            spike_out_1 : out std_logic;
            voltage_out_0 : out integer;
            voltage_out_1 : out integer
        );
    end component;

    -- Signals for testbench
    signal clk_tb         : std_logic := '0';
    signal reset_n_tb     : std_logic := '0';
    signal sp_0_tb        : std_logic := '0';
    signal sp_1_tb        : std_logic := '0';
    signal sp_2_tb        : std_logic := '0';
    signal sp_3_tb        : std_logic := '0';
    signal sp_4_tb        : std_logic := '0';
    signal sp_5_tb        : std_logic := '0';
    signal sp_6_tb        : std_logic := '0';
    signal sp_7_tb        : std_logic := '0';
    signal sp_8_tb        : std_logic := '0';
    signal sp_9_tb        : std_logic := '0';
    signal spike_out_0_tb : std_logic;
    signal spike_out_1_tb : std_logic;
    signal voltage_out_0_tb : integer;
    signal voltage_out_1_tb : integer;

begin
    -- Instantiate the Unit Under Test (UUT)
    uut: snn_spikes
        port map (
            clk => clk_tb,
            reset_n => reset_n_tb,
            sp_0 => sp_0_tb,
            sp_1 => sp_1_tb,
            sp_2 => sp_2_tb,
            sp_3 => sp_3_tb,
            sp_4 => sp_4_tb,
            sp_5 => sp_5_tb,
            sp_6 => sp_6_tb,
            sp_7 => sp_7_tb,
            sp_8 => sp_8_tb,
            sp_9 => sp_9_tb,
            spike_out_0 => spike_out_0_tb,
            spike_out_1 => spike_out_1_tb,
            voltage_out_0 => voltage_out_0_tb,
            voltage_out_1 => voltage_out_1_tb
        );

    -- Clock generation process
    clk_process : process
    begin
        clk_tb <= '0';
        wait for 10 ns;
        clk_tb <= '1';
        wait for 10 ns;
    end process;

    -- Stimuli process
    stim_proc: process
    begin        
        -- Apply reset
        reset_n_tb <= '0';
        wait for 20 ns;
        reset_n_tb <= '1';
        wait for 20 ns;
        
        for i in 1 to 3 loop
            -- Test case 1: Apply some spike signals and enable input
            sp_0_tb <= '1';
            sp_1_tb <= '0';
            sp_2_tb <= '1';
            sp_3_tb <= '0';
            sp_4_tb <= '1';
            sp_5_tb <= '0';
            sp_6_tb <= '1';
            sp_7_tb <= '0';
            sp_8_tb <= '1';
            sp_9_tb <= '0';
            wait for 10 ns;
    
            -- Test case 2: Change spike signals
            sp_0_tb <= '0';
            sp_1_tb <= '1';
            sp_2_tb <= '0';
            sp_3_tb <= '1';
            sp_4_tb <= '0';
            sp_5_tb <= '1';
            sp_6_tb <= '0';
            sp_7_tb <= '1';
            sp_8_tb <= '0';
            sp_9_tb <= '1';
            wait for 10 ns;
    
            -- Test case 3: No spikes and reset
            sp_0_tb <= '0';
            sp_1_tb <= '0';
            sp_2_tb <= '0';
            sp_3_tb <= '0';
            sp_4_tb <= '0';
            sp_5_tb <= '0';
            sp_6_tb <= '0';
            sp_7_tb <= '0';
            sp_8_tb <= '0';
            sp_9_tb <= '0';
            wait for 10 ns;
    
            -- Test case 4: Additional pattern and testing
            sp_0_tb <= '1';
            sp_1_tb <= '1';
            sp_2_tb <= '0';
            sp_3_tb <= '0';
            sp_4_tb <= '0';
            sp_5_tb <= '1';
            sp_6_tb <= '1';
            sp_7_tb <= '0';
            sp_8_tb <= '1';
            sp_9_tb <= '0';
            wait for 10 ns;
        end loop;
        -- End of simulation
        wait for 10 ns;
    end process;

end behavior;