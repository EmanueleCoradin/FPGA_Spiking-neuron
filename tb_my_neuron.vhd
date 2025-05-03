library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity tb_neuron is
end tb_neuron;

architecture behavior of tb_neuron is

    -- Signals
    signal clk         : std_logic := '0';
    signal reset       : std_logic := '1';

    signal sp_0        : std_logic := '0';
    signal sp_1        : std_logic := '0';
    signal sp_2        : std_logic := '0';
    signal sp_3        : std_logic := '0';
    signal sp_4        : std_logic := '0';
    signal sp_5        : std_logic := '0';
    signal sp_6        : std_logic := '0';
    signal sp_7        : std_logic := '0';
    signal sp_8        : std_logic := '0';
    signal sp_9        : std_logic := '0';

    signal spike_out   : std_logic;
    signal voltage_out : integer;

    constant clk_period : time := 10 ns;

begin

    -- DUT instantiation
    uut: entity work.neuron
        generic map(
            w_0 => 10, w_1 => 5, w_2 => 3, w_3 => 8, w_4 => 1,
            w_5 => 10, w_6 => 7, w_7 => 3, w_8 => 5, w_9 => 4,
            v_th => 5,
            v_baseline => -7,
            decay_factor => 9,
            decay_divisor => 10,
            N_0 => 2, N_1 => 4, N_2 => 1, N_3 => 7, N_4 => 5,
            N_5 => 9, N_6 => 6, N_7 => 2, N_8 => 4, N_9 => 8
        )
        port map (
            clk => clk,
            reset => reset,
            sp_0 => sp_0, sp_1 => sp_1, sp_2 => sp_2,
            sp_3 => sp_3, sp_4 => sp_4, sp_5 => sp_5,
            sp_6 => sp_6, sp_7 => sp_7, sp_8 => sp_8,
            sp_9 => sp_9,
            spike_out => spike_out,
            voltage_out => voltage_out
        );

    -- Clock generation
    clk_process : process
    begin
        while true loop
            clk <= '0';
            wait for clk_period / 2;
            clk <= '1';
            wait for clk_period / 2;
        end loop;
    end process;

    -- Logging on each clock edge
    log_process : process(clk)
    begin
        if rising_edge(clk) then
            report "Voltage: " & integer'image(voltage_out) &
                   " | Spike: " & std_logic'image(spike_out);
        end if;
    end process;

    -- Stimulus
    stim_proc: process
    begin
        -- Initial reset
        wait for 20 ns;
        reset <= '0';

        -- Stimuli pattern
        wait for 20 ns;
        sp_0 <= '1'; wait for clk_period;
        sp_0 <= '0';

        wait for 30 ns;
        sp_1 <= '1'; wait for clk_period;
        sp_1 <= '0';

        wait for 20 ns;
        sp_2 <= '1'; wait for clk_period;
        sp_2 <= '0';

        wait for 200 ns;
        sp_0 <= '1'; wait for clk_period;
        sp_0 <= '0';

        -- Extend simulation time
        wait for 400 ns;

        -- Done
        report "Testbench finished.";
        wait;
    end process;

end behavior;
