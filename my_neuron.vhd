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
		w_0		      : integer:= 0;
		w_1		      : integer:= 0;
		w_2		      : integer:= 0;
		w_3		      : integer:= 0;
		w_4		      : integer:= 0;
		w_5		      : integer:= 0;
		w_6		      : integer:= 0;
		w_7		      : integer:= 0;
		w_8		      : integer:= 0;
		w_9		      : integer:= 0;
		
		
		v_th	      : integer:= 0;    -- Membrane threshold voltage
		v_baseline    : integer:= 0;    -- Negative value to which the potential is reset once it reaches the threshold
		decay_factor  : integer := 9;   -- Decay factor for leakage
		decay_divisor : integer := 10;  -- Decay divisor for leakage
		
		-- Length of the shift registers for the delays -- INITALIZED TO 64 or 32
		N_0     : integer:= 0;
		N_1     : integer:= 0;
		N_2     : integer:= 0;
		N_3     : integer:= 0;
		N_4     : integer:= 0;
		N_5     : integer:= 0;
		N_6     : integer:= 0;
		N_7     : integer:= 0;
		N_8     : integer:= 0;
		N_9     : integer:= 0
	);

	port(
		clk			: in  std_logic := '0';  -- Clock input
		reset		: in  std_logic := '0';  -- Global reset
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
		spike_out	: out std_logic := '0';   -- Output spike signal
		voltage_out : out integer   :=  0
	);
end neuron;

architecture behave of neuron is

    -- If these internal signals reach 1, then add the corresponding weight
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

	-- Signals for adder tree to sum inputs more efficiently
	signal tmp_sum_0_1	           : integer:= 0;
	signal tmp_sum_2_3	           : integer:= 0;
	signal tmp_sum_4_5	           : integer:= 0;
	signal tmp_sum_6_7	           : integer:= 0;
	signal tmp_sum_8_9	           : integer:= 0;
	signal tmp_sum_0_1_2_3         : integer:= 0;
	signal tmp_sum_4_5_6_7         : integer:= 0;
	signal tmp_sum_0_1_2_3_4_5_6_7 : integer:= 0;
	signal sum	                   : integer:= 0;
	
	-- Memory for the spike; 
	-- Useful for registering the voltage throughout the time
	signal is_spike : boolean := FALSE;

begin

    delay0 : entity work.shift_register
        Generic map(N => N_0)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_0,
                    output => out_0);
                    
    delay1 : entity work.shift_register
        Generic map(N => N_1)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_1,
                    output => out_1);
                    
    delay2 : entity work.shift_register
        Generic map(N => N_2)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_2,
                    output => out_2);
                    
    delay3 : entity work.shift_register
        Generic map(N => N_3)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_3,
                    output => out_3);
                   
    delay4 : entity work.shift_register
        Generic map(N => N_4)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_4,
                    output => out_4);
                   
    delay5 : entity work.shift_register
        Generic map(N => N_5)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_5,
                    output => out_5);
                    
    delay6 : entity work.shift_register
        Generic map(N => N_6)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_6,
                    output => out_6);
                    
    delay7 : entity work.shift_register
        Generic map(N => N_7)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_7,
                    output => out_7);
                    
    delay8 : entity work.shift_register
        Generic map(N => N_8)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_8,
                    output => out_8);
                    
    delay9 : entity work.shift_register
        Generic map(N => N_9)
        Port map   (Clk => clk,
                    Rst => reset,
                    input => sp_9,
                    output => out_9);

process
	-- Internal voltage of the neuron (membrane potential)
	variable voltage : integer := 0;	

begin
	wait until rising_edge(clk);  -- Wait for clock edge
	
	-- Reset if the flag is active
	if reset = '1' then
	   voltage := 0;
	   spike_out <= '0';
	   is_spike <= FALSE;

    else
        -- Assign weight to tmp_sum_X if corresponding spike is high, else 0
        if out_0 = '1' then
            tmp_sum_0 <= w_0;
        else
            tmp_sum_0 <= 0;
        end if;
    
        if out_1 = '1' then
            tmp_sum_1 <= w_1;
        else
            tmp_sum_1 <= 0;
        end if;
    
        if out_2 = '1' then
            tmp_sum_2 <= w_2;
        else
            tmp_sum_2 <= 0;
        end if;
    
        if out_3 = '1' then
            tmp_sum_3 <= w_3;
        else
            tmp_sum_3 <= 0;
        end if;
    
        if out_4 = '1' then
            tmp_sum_4 <= w_4;
        else
            tmp_sum_4 <= 0;
        end if;
    
        if out_5 = '1' then
            tmp_sum_5 <= w_5;
        else
            tmp_sum_5 <= 0;
        end if;
    
        if out_6 = '1' then
            tmp_sum_6 <= w_6;
        else
            tmp_sum_6 <= 0;
        end if;
        
        if out_7 = '1' then
            tmp_sum_7 <= w_7;
        else
            tmp_sum_7 <= 0;
        end if;
        
        if out_8 = '1' then
            tmp_sum_8 <= w_8;
        else
            tmp_sum_8 <= 0;
        end if;
        
        if out_9 = '1' then
            tmp_sum_9 <= w_9;
        else
            tmp_sum_9 <= 0;
        end if;
    
        -- Efficient adder tree to compute final weighted sum
        tmp_sum_0_1             <= tmp_sum_0 + tmp_sum_1;
        tmp_sum_2_3             <= tmp_sum_2 + tmp_sum_3;
        tmp_sum_4_5             <= tmp_sum_4 + tmp_sum_5;
        tmp_sum_6_7             <= tmp_sum_6 + tmp_sum_7;
        tmp_sum_8_9             <= tmp_sum_8 + tmp_sum_9;
        
        tmp_sum_0_1_2_3         <= tmp_sum_0_1 + tmp_sum_2_3;
        tmp_sum_4_5_6_7         <= tmp_sum_4_5 + tmp_sum_6_7;
        
        tmp_sum_0_1_2_3_4_5_6_7 <= tmp_sum_0_1_2_3 + tmp_sum_4_5_6_7;
        
        sum                     <= tmp_sum_0_1_2_3_4_5_6_7 + tmp_sum_8_9;
    
        -- Integrate weighted input sum into membrane potential with leakage
        voltage := (voltage * decay_factor) / decay_divisor + sum;
    
        if is_spike then
        -- Spike in previous cycle -> Set the voltage to v_baseline <= 0
            voltage := v_baseline;
            is_spike <= FALSE;

        elsif voltage > v_th then
        -- Check if voltage exceeds threshold -> fire a spike
            spike_out <= '1';
            is_spike <= TRUE;               
        else
        -- No spike in this or previous cycle -> Keep working
            spike_out <= '0';           -- No spike
        end if;
    end if;
    
	
    voltage_out <= voltage;

end process;

end behave;
