library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity snn_architecture is
    port (
        clk         : in  std_logic;                         -- Input clock
        reset_n     : in  std_logic;                         -- Neuron reset

        -- Spike train inputs (one-hot encoded spike signals)
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

        -- Output signals from neuron 0
        spike_out_0    : out std_logic;                      -- Output spike
        voltage_out_0  : out std_logic_vector(7 downto 0)    -- Membrane voltage
    );
end snn_architecture;

architecture behave of snn_architecture is

    -----------------------------------------------------------------------
    -- Internal signals
    -----------------------------------------------------------------------
    signal output_spike_0 : std_logic;                        -- Neuron spike output
    signal v_out_0        : std_logic_vector(7 downto 0);     -- Neuron membrane voltage

    -----------------------------------------------------------------------
    -- Neuron parameters
    -----------------------------------------------------------------------
    constant v_th           : integer := 100;                 -- Threshold voltage
    constant v_baseline     : integer := -30;                 -- Baseline/resting potential
    constant decay_factor   : integer := 90;                  -- Decay numerator
    constant decay_divisor  : integer := 100;                 -- Decay denominator

begin

    -----------------------------------------------------------------------
    -- Instantiation of output neuron 0
    -----------------------------------------------------------------------
    output_neuron_0 : entity work.neuron
        generic map(
            -- Synaptic weights for each input feature
            w_0 => 13, w_1 => 17, w_2 => 19, w_3 => 18, w_4 => 15,
            w_5 => 18, w_6 => 14, w_7 => 19, w_8 => 14, w_9 => 15,

            -- Voltage and decay parameters
            v_th         => v_th,
            v_baseline   => v_baseline,
            decay_factor => decay_factor,
            decay_divisor => decay_divisor,

            -- Neuron delays
            N_0 => 80, N_1 => 77, N_2 => 74, N_3 => 70, N_4 => 62,
            N_5 => 53, N_6 => 41, N_7 => 26, N_8 => 10, N_9 => 1
        )
        port map(
            clk         => clk,
            reset       => reset_n,
            sp_0        => sp_0, sp_1 => sp_1, sp_2 => sp_2, sp_3 => sp_3, sp_4 => sp_4,
            sp_5        => sp_5, sp_6 => sp_6, sp_7 => sp_7, sp_8 => sp_8, sp_9 => sp_9,
            spike_out   => output_spike_0,
            voltage_out => v_out_0
        );

    -----------------------------------------------------------------------
    -- Drive output ports from internal neuron signals
    -----------------------------------------------------------------------
    spike_out_0   <= output_spike_0;
    voltage_out_0 <= v_out_0;

end behave;
