-- Testbench for the neuron VHDL model
-- Author: Your Name
-- This testbench stimulates the neuron entity and checks its behavior

library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity tb_neuron is
    -- No ports in a testbench
end tb_neuron;

architecture behavior of tb_neuron is

    -- Signal declarations to connect to the neuron entity
    signal clk            : std_logic := '0';
    signal reset          : std_logic := '0';
    signal sp_0           : std_logic := '0';
    signal sp_1           : std_logic := '0';
    signal sp_2           : std_logic := '0';
    signal sp_3           : std_logic := '0';
    signal sp_4           : std_logic := '0';
    signal sp_5           : std_logic := '0';
    signal sp_6           : std_logic := '0';
    signal neuron_reset   : std_logic := '0';
    signal spike_out      : std_logic := '0';
    signal voltage_out    : integer   :=  0 ;

    -- Constants for neuron weights and bias
    constant w_0 : integer := 10;
    constant w_1 : integer := 20;
    constant w_2 : integer := 30;
    constant w_3 : integer := 40;
    constant w_4 : integer := 50;
    constant w_5 : integer := 60;
    constant w_6 : integer := 70;
    constant bias : integer := 5;
    constant v_th : integer := 100;

begin
    -- Instantiate the neuron entity
    uut: entity work.neuron
        generic map (
            w_0   => w_0,
            w_1   => w_1,
            w_2   => w_2,
            w_3   => w_3,
            w_4   => w_4,
            w_5   => w_5,
            w_6   => w_6,
            bias  => bias,
            v_th  => v_th
        )
        port map (
            clk         => clk,
            reset       => reset,
            sp_0        => sp_0,
            sp_1        => sp_1,
            sp_2        => sp_2,
            sp_3        => sp_3,
            sp_4        => sp_4,
            sp_5        => sp_5,
            sp_6        => sp_6,
            neuron_reset => neuron_reset,
            spike_out   => spike_out,
            voltage_out => voltage_out
        );

    -- Generate clock signal
    clk_process :process
    begin
        clk <= not clk;
        wait for 10 ns;
    end process;

    -- Stimulus process
    stim_proc: process
    begin
        -- Initialize signals
        reset <= '0';
        sp_0 <= '0';
        sp_1 <= '0';
        sp_2 <= '0';
        sp_3 <= '0';
        sp_4 <= '0';
        sp_5 <= '0';
        sp_6 <= '0';
        neuron_reset <= '0';
        wait for 20 ns;
        
        for i in 0 to 3 loop
            
            -- Reset neuron and apply spikes
            sp_0 <= '1'; -- simulate spike at input 0
            wait for 20 ns;
            
            sp_0 <= '0'; -- stop spike at input 0
            sp_1 <= '1'; -- simulate spike at input 1
            wait for 20 ns;
            
            sp_1 <= '0'; -- stop spike at input 1
            sp_2 <= '1'; -- simulate spike at input 2
            wait for 20 ns;
            
            -- Apply additional spikes if necessary
            sp_2 <= '0'; 
            sp_3 <= '1';
            wait for 20 ns;
            
            sp_3 <= '0'; 
            sp_4 <= '1';
            wait for 20 ns;
            
            sp_4 <= '0';
            wait for 20 ns;
            
            sp_0 <= '1';
            sp_1 <= '1';
            wait for 20 ns;
            
            sp_0 <= '0';
            sp_1 <= '0';
            wait for 20 ns;
        end loop;
        
        -- End simulation after a certain amount of time
        wait for 100 ns;
        assert false report "Simulation Finished" severity failure;
    end process;

end behavior;