-- neuron.vhd
-- Neuron for Spiking Neural Network
--
-- Author: Klaus Niederberger
-- Release: Marco Winzker, Hochschule Bonn-Rhein-Sieg, 22.12.2022
-- FPGA Vision Remote Lab http://h-brs.de/fpga-vision-lab

library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

-- Entity declaration of a spiking neuron
-- Takes spike inputs and weights, computes weighted sum,
-- integrates over time, and generates an output spike if threshold is exceeded.
entity neuron is
	generic(
		-- Input weights for each synapse
		w_0		: integer:= 0;
		w_1		: integer:= 0;
		w_2		: integer:= 0;
		w_3		: integer:= 0;
		w_4		: integer:= 0;
		w_5		: integer:= 0;
		w_6		: integer:= 0;
		w_7		: integer:= 0;
		w_8		: integer:= 0;
		w_9		: integer:= 0;
		w_10    : integer:= 0;
		w_11    : integer:= 0;
		w_12    : integer:= 0;
		w_13    : integer:= 0;
		w_14    : integer:= 0;
		bias	: integer:= 0;   -- Bias added to the input sum
		v_th	: integer:= 0    -- Membrane threshold voltage
		-- decay_factor : integer := 9;  -- Decay factor for leakage
		-- decay_divisor : integer := 10  -- Decay divisor for leakage
	);

	port(
		clk			: in  std_logic := '0';  -- Clock input
		reset		: in  std_logic := '0';  -- Global reset (unused here)
		sp_0		: in  std_logic := '0';  -- Spike inputs (1-bit each)
		sp_1		: in  std_logic := '0';
		sp_2		: in  std_logic := '0';
		sp_3		: in  std_logic := '0';
		sp_4		: in  std_logic := '0';
		sp_5		: in  std_logic := '0';
		sp_6		: in  std_logic := '0';
		sp_7		: in  std_logic := '0';
		sp_8		: in  std_logic := '0';
		sp_9		: in  std_logic := '0';
		sp_10		: in  std_logic := '0';
		sp_11		: in  std_logic := '0';
		sp_12		: in  std_logic := '0';
		sp_13		: in  std_logic := '0';
		sp_14		: in  std_logic := '0';
		neuron_reset: in  std_logic := '0';  -- Resets internal voltage to 0
		spike_out	: out std_logic := '0';   -- Output spike signal
		voltage_out : out integer   :=   0
	);
end neuron;

architecture behave of neuron is

	-- Intermediate signals for conditional weight accumulation
	signal tmp_sum_0			: integer:= 0;
	signal tmp_sum_1			: integer:= 0;
	signal tmp_sum_2			: integer:= 0;
	signal tmp_sum_3			: integer:= 0;
	signal tmp_sum_4			: integer:= 0;
	signal tmp_sum_5			: integer:= 0;
	signal tmp_sum_6			: integer:= 0;
	signal tmp_sum_7			: integer:= 0;
	signal tmp_sum_8			: integer:= 0;
	signal tmp_sum_9			: integer:= 0;
	signal tmp_sum_10			: integer:= 0;
    signal tmp_sum_11			: integer:= 0;
    signal tmp_sum_12			: integer:= 0;
    signal tmp_sum_13			: integer:= 0;
    signal tmp_sum_14			: integer:= 0;


	-- Signals for adder tree to sum inputs more efficiently
	signal tmp_sum_b_0	       : integer:= 0;
	signal tmp_sum_1_2	       : integer:= 0;
	signal tmp_sum_3_4	       : integer:= 0;
	signal tmp_sum_5_6	       : integer:= 0;
	signal tmp_sum_7_8	       : integer:= 0;
	signal tmp_sum_9_10	       : integer:= 0;
	signal tmp_sum_11_12       : integer:= 0;
	signal tmp_sum_13_14       : integer:= 0;
	signal tmp_sum_b_0_1_2	   : integer:= 0;
	signal tmp_sum_3_4_5_6	   : integer:= 0;
	signal tmp_sum_7_8_9_10	   : integer:= 0;
	signal tmp_sum_11_12_13_14 : integer:= 0;
	signal tmp_sum_half_1      : integer:= 0;
	signal tmp_sum_half_2      : integer:= 0;
	signal sum	               : integer:= 0;

begin

process
	-- Internal voltage of the neuron (membrane potential)
	variable voltage : integer := 0;	

begin
	wait until rising_edge(clk);  -- Wait for clock edge

	-- Assign weight to tmp_sum_X if corresponding spike is high, else 0
	if sp_0 = '1' then
		tmp_sum_0 <= w_0;
	else
		tmp_sum_0 <= 0;
	end if;

	if sp_1 = '1' then
		tmp_sum_1 <= w_1;
	else
		tmp_sum_1 <= 0;
	end if;

	if sp_2 = '1' then
		tmp_sum_2 <= w_2;
	else
		tmp_sum_2 <= 0;
	end if;

	if sp_3 = '1' then
		tmp_sum_3 <= w_3;
	else
		tmp_sum_3 <= 0;
	end if;

	if sp_4 = '1' then
		tmp_sum_4 <= w_4;
	else
		tmp_sum_4 <= 0;
	end if;

	if sp_5 = '1' then
		tmp_sum_5 <= w_5;
	else
		tmp_sum_5 <= 0;
	end if;

	if sp_6 = '1' then
		tmp_sum_6 <= w_6;
	else
		tmp_sum_6 <= 0;
	end if;
	
	if sp_7 = '1' then
		tmp_sum_7 <= w_7;
	else
		tmp_sum_7 <= 0;
	end if;
	
	if sp_8 = '1' then
		tmp_sum_8 <= w_8;
	else
		tmp_sum_8 <= 0;
	end if;
	
	if sp_9 = '1' then
		tmp_sum_9 <= w_9;
	else
		tmp_sum_9 <= 0;
	end if;
	
	if sp_10 = '1' then
		tmp_sum_10 <= w_10;
	else
		tmp_sum_10 <= 0;
	end if;
	
    if sp_11 = '1' then
		tmp_sum_11 <= w_11;
	else
		tmp_sum_11 <= 0;
	end if;
	
    if sp_12 = '1' then
		tmp_sum_12 <= w_12;
	else
		tmp_sum_12 <= 0;
	end if;
	
    if sp_13 = '1' then
		tmp_sum_13 <= w_13;
	else
        tmp_sum_13 <= 0;
	end if;
	
    if sp_14 = '1' then
		tmp_sum_14 <= w_14;
	else
		tmp_sum_14 <= 0;
	end if;

	-- Efficient adder tree to compute final weighted sum
	tmp_sum_b_0      <= bias        + tmp_sum_0;
	tmp_sum_1_2      <= tmp_sum_1   + tmp_sum_2;
	tmp_sum_3_4      <= tmp_sum_3   + tmp_sum_4;
	tmp_sum_5_6      <= tmp_sum_5   + tmp_sum_6;
	tmp_sum_7_8      <= tmp_sum_7   + tmp_sum_8;
	tmp_sum_9_10     <= tmp_sum_9   + tmp_sum_10;
	tmp_sum_11_12    <= tmp_sum_11  + tmp_sum_12;
	tmp_sum_13_14    <= tmp_sum_13  + tmp_sum_14;

	tmp_sum_b_0_1_2     <= tmp_sum_b_0 + tmp_sum_1_2;
	tmp_sum_3_4_5_6     <= tmp_sum_3_4 + tmp_sum_5_6;
	tmp_sum_7_8_9_10    <= tmp_sum_7_8 + tmp_sum_9_10;
	tmp_sum_11_12_13_14 <= tmp_sum_11_12 + tmp_sum_13_14;

	tmp_sum_half_1 <= tmp_sum_b_0_1_2 + tmp_sum_3_4_5_6;
	tmp_sum_half_2 <= tmp_sum_7_8_9_10 + tmp_sum_11_12_13_14;
	
	sum <= tmp_sum_half_1 + tmp_sum_half_2;

	-- Integrate weighted input sum into membrane potential
	voltage := voltage + sum;

	-- Integrate weighted input sum into membrane potential with leakage
	-- voltage := (voltage * decay_factor) / decay_divisor + sum;

	-- Check if voltage exceeds threshold -> fire a spike
	if voltage > v_th then
		voltage := voltage - v_th;  -- Reset voltage by subtracting threshold
		spike_out <= '1';
	else
		spike_out <= '0';           -- No spike
	end if;

	-- External reset input resets membrane voltage
	if neuron_reset = '1' then
		voltage := 0;
	end if;
	
    voltage_out <= voltage;

end process;

end behave;
