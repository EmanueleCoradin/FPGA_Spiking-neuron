library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

-- Entity declaration of a spiking neuron
-- Processes spike inputs with weights, integrates over time,
-- applies leakage decay, and generates output spike if threshold crossed.
entity neuron is
    generic (
        -- Input synaptic weights for each input spike
        w_0           : integer := 0;
        w_1           : integer := 0;
        w_2           : integer := 0;
        w_3           : integer := 0;
        w_4           : integer := 0;
        w_5           : integer := 0;
        w_6           : integer := 0;
        w_7           : integer := 0;
        w_8           : integer := 0;
        w_9           : integer := 0;

        -- Neuron parameters
        v_th          : integer := 100;    -- Threshold voltage for firing
        v_baseline    : integer := -20;    -- Reset voltage after spike
        decay_factor  : integer := 39;     -- Numerator for voltage decay factor
        decay_divisor : integer := 40;     -- Denominator for voltage decay factor

        -- Delay lengths for each input's shift register (for synaptic delay modeling)
        N_0           : integer := 0;
        N_1           : integer := 0;
        N_2           : integer := 0;
        N_3           : integer := 0;
        N_4           : integer := 0;
        N_5           : integer := 0;
        N_6           : integer := 0;
        N_7           : integer := 0;
        N_8           : integer := 0;
        N_9           : integer := 0
    );

    port (
        clk          : in  std_logic := '0';    -- Clock input
        reset        : in  std_logic := '0';    -- Asynchronous reset (active high)

        -- Spike inputs (one bit per input channel)
        sp_0         : in  std_logic := '0';
        sp_1         : in  std_logic := '0';
        sp_2         : in  std_logic := '0';
        sp_3         : in  std_logic := '0';
        sp_4         : in  std_logic := '0';
        sp_5         : in  std_logic := '0';
        sp_6         : in  std_logic := '0';
        sp_7         : in  std_logic := '0';
        sp_8         : in  std_logic := '0';
        sp_9         : in  std_logic := '0';

        spike_out    : out std_logic := '0';    -- Output spike signal
        voltage_out  : out std_logic_vector(7 downto 0)  -- Membrane voltage output (8-bit signed)
    );
end neuron;

architecture behave of neuron is

    -- Lookup table for voltage decay based on current voltage value
    type int_array is array (0 to 511) of integer;
    signal decay_LUT : int_array;

    -- Internal delayed versions of spike inputs (shift registers)
    signal out_0 : std_logic := '0';
    signal out_1 : std_logic := '0';
    signal out_2 : std_logic := '0';
    signal out_3 : std_logic := '0';
    signal out_4 : std_logic := '0';
    signal out_5 : std_logic := '0';
    signal out_6 : std_logic := '0';
    signal out_7 : std_logic := '0';
    signal out_8 : std_logic := '0';
    signal out_9 : std_logic := '0';

    -- Temporary sums for weights conditional on spikes (initial stage)
    signal tmp_sum_0 : integer := 0;
    signal tmp_sum_1 : integer := 0;
    signal tmp_sum_2 : integer := 0;
    signal tmp_sum_3 : integer := 0;
    signal tmp_sum_4 : integer := 0;
    signal tmp_sum_5 : integer := 0;
    signal tmp_sum_6 : integer := 0;
    signal tmp_sum_7 : integer := 0;
    signal tmp_sum_8 : integer := 0;
    signal tmp_sum_9 : integer := 0;

    -- Additional signals for hierarchical adder tree
    signal tmp_sum_a  : integer := 0;
    signal tmp_sum_b  : integer := 0;
    signal tmp_sum_c  : integer := 0;
    signal tmp_sum_d  : integer := 0;
    signal tmp_sum_e  : integer := 0;
    signal tmp_sum_f  : integer := 0;

    -- Intermediate sums for adder tree stages
    signal tmp_sum_0_1             : integer := 0;
    signal tmp_sum_2_3             : integer := 0;
    signal tmp_sum_4_5             : integer := 0;
    signal tmp_sum_6_7             : integer := 0;
    signal tmp_sum_8_9             : integer := 0;
    signal tmp_sum_a_b             : integer := 0;
    signal tmp_sum_c_d             : integer := 0;
    signal tmp_sum_e_f             : integer := 0;
    signal tmp_sum_0_1_2_3         : integer := 0;
    signal tmp_sum_4_5_6_7         : integer := 0;
    signal tmp_sum_8_9_a_b         : integer := 0;
    signal tmp_sum_c_d_e_f         : integer := 0;
    signal tmp_sum_0_1_2_3_4_5_6_7 : integer := 0;
    signal tmp_sum_8_9_a_b_c_d_e_f : integer := 0;

    -- Final weighted sum of all inputs
    signal sum : integer := 0;

    -- Membrane potential and decayed membrane potential
    signal voltage         : integer := 0;
    signal decayed_voltage : integer := 0;

begin

    ---------------------------------------------------------------------------
    -- Instantiate delay shift registers for each input spike signal
    -- Models synaptic delay by delaying input spikes by N_i cycles
    ---------------------------------------------------------------------------
    delay0 : entity work.shift_register
        generic map (N => N_0)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_0,
            output=> out_0
        );

    delay1 : entity work.shift_register
        generic map (N => N_1)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_1,
            output=> out_1
        );

    delay2 : entity work.shift_register
        generic map (N => N_2)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_2,
            output=> out_2
        );

    delay3 : entity work.shift_register
        generic map (N => N_3)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_3,
            output=> out_3
        );

    delay4 : entity work.shift_register
        generic map (N => N_4)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_4,
            output=> out_4
        );

    delay5 : entity work.shift_register
        generic map (N => N_5)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_5,
            output=> out_5
        );

    delay6 : entity work.shift_register
        generic map (N => N_6)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_6,
            output=> out_6
        );

    delay7 : entity work.shift_register
        generic map (N => N_7)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_7,
            output=> out_7
        );

    delay8 : entity work.shift_register
        generic map (N => N_8)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_8,
            output=> out_8
        );

    delay9 : entity work.shift_register
        generic map (N => N_9)
        port map (
            Clk   => clk,
            Rst   => reset,
            input => sp_9,
            output=> out_9
        );

    ---------------------------------------------------------------------------
    -- Main process: handles weighted sum, integration, decay, and spike generation
    ---------------------------------------------------------------------------
    process(clk)
    begin
        -- Fill decay lookup table 
        -- This LUT maps voltage values to decayed values, modeling leakage
        for iv in 0 to 128 + v_th loop
            decay_LUT(iv) <= ((iv - 128) * decay_factor) / decay_divisor;
        end loop;

        for iv in 128 + v_th + 1 to 511 loop
            decay_LUT(iv) <= 0;
        end loop;

        -- Apply voltage decay lookup
        decayed_voltage <= decay_LUT(voltage + 128);
        
        if rising_edge(clk) then

            -- Handle asynchronous reset
            if reset = '1' then
                voltage   <= 0;
                spike_out <= '0';

            else
                -- Assign weights based on delayed spike inputs (if spike present, assign weight, else zero)
                if out_0 = '1' then tmp_sum_0 <= w_0; else tmp_sum_0 <= 0; end if;
                if out_1 = '1' then tmp_sum_1 <= w_1; else tmp_sum_1 <= 0; end if;
                if out_2 = '1' then tmp_sum_2 <= w_2; else tmp_sum_2 <= 0; end if;
                if out_3 = '1' then tmp_sum_3 <= w_3; else tmp_sum_3 <= 0; end if;
                if out_4 = '1' then tmp_sum_4 <= w_4; else tmp_sum_4 <= 0; end if;
                if out_5 = '1' then tmp_sum_5 <= w_5; else tmp_sum_5 <= 0; end if;
                if out_6 = '1' then tmp_sum_6 <= w_6; else tmp_sum_6 <= 0; end if;
                if out_7 = '1' then tmp_sum_7 <= w_7; else tmp_sum_7 <= 0; end if;
                if out_8 = '1' then tmp_sum_8 <= w_8; else tmp_sum_8 <= 0; end if;
                if out_9 = '1' then tmp_sum_9 <= w_9; else tmp_sum_9 <= 0; end if;
                
                -- Combine sums using hierarchical adder tree for efficiency
                tmp_sum_0_1             <= tmp_sum_0 + tmp_sum_1;
                tmp_sum_2_3             <= tmp_sum_2 + tmp_sum_3;
                tmp_sum_4_5             <= tmp_sum_4 + tmp_sum_5;
                tmp_sum_6_7             <= tmp_sum_6 + tmp_sum_7;
                tmp_sum_8_9             <= tmp_sum_8 + tmp_sum_9;

                -- Intermediate sums a,b,c,d,e,f are zero by default; can be extended if needed
                tmp_sum_a_b             <= tmp_sum_a + tmp_sum_b;
                tmp_sum_c_d             <= tmp_sum_c + tmp_sum_d;
                tmp_sum_e_f             <= tmp_sum_e + tmp_sum_f;

                tmp_sum_0_1_2_3         <= tmp_sum_0_1 + tmp_sum_2_3;
                tmp_sum_4_5_6_7         <= tmp_sum_4_5 + tmp_sum_6_7;
                tmp_sum_8_9_a_b         <= tmp_sum_8_9 + tmp_sum_a_b;
                tmp_sum_c_d_e_f         <= tmp_sum_c_d + tmp_sum_e_f;

                tmp_sum_0_1_2_3_4_5_6_7 <= tmp_sum_0_1_2_3 + tmp_sum_4_5_6_7;
                tmp_sum_8_9_a_b_c_d_e_f <= tmp_sum_8_9_a_b + tmp_sum_c_d_e_f;

                -- Final sum of all weighted spikes
                sum <= tmp_sum_0_1_2_3_4_5_6_7 + tmp_sum_8_9_a_b_c_d_e_f;

                -- Integrate input (weighted sum) with decayed voltage
                voltage <= decayed_voltage + sum;

                -- Spike generation and reset logic
                if voltage >= v_th then
                    voltage   <= v_baseline;  -- Reset voltage to baseline after spike
                    spike_out <= '1';         -- Output spike generated
                else
                    spike_out <= '0';         -- No spike output
                end if;

            end if;

            -- Clamp voltage output to signed 8-bit range (-128 to 127)
            if voltage > 127 then
                voltage_out <= std_logic_vector(to_signed(127, 8));
            else
                voltage_out <= std_logic_vector(to_signed(voltage, 8));
            end if;

        end if;
    end process;

end behave;
