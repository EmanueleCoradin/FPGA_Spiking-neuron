library ieee;
use ieee.std_logic_1164.all;  -- Standard logic library for digital design
use ieee.numeric_std.all;     -- Numeric library for unsigned and signed types

-- Entity declaration: Defines the interface of the sampler_generator component
entity sampler_generator is
  port (
    clock        : in  std_logic;         -- Input clock signal
    uart_rx      : in  std_logic;         -- UART receive data signal (to be sampled)
    baudrate_out : out std_logic          -- Output signal indicating the baud rate timing
  );
end entity sampler_generator;

-- Architecture: Describes the internal implementation of the component (RTL - Register Transfer Level)
architecture rtl of sampler_generator is

  -- Define the states for the state machine (Idle, Start, Bit1 to Bit8, Stop)
  type state_t is (idle_s, start_s, bit0_s, bit1_s, bit2_s, bit3_s, bit4_s, bit5_s, bit6_s, bit7_s, bit8_s, stop_s);
  signal state : state_t := idle_s;  -- Initial state is idle

  -- Signal definitions for counters and control
  signal counter        : unsigned(10 downto 0) := (others => '0');      -- 11-bit counter for baud rate generation
  signal delay_counter  : unsigned(10 downto 0) := (others => '0');      -- 11-bit delay counter for baud rate
  constant divisor      : unsigned(10 downto 0) := to_unsigned(867, 11); -- Divisor to control baud rate (for example, 867 clock cycles)
  constant half_divisor : unsigned(10 downto 0) := to_unsigned(433, 11); -- Half divisor for generating correct baud rate timing
  signal busy           : std_logic             := '0';                   -- Indicates if the system is busy
  signal pulse_out      : std_logic;                                            -- Output pulse for baud rate
  signal enable_counter : std_logic             := '0';                   -- Enable signal for the counter
  signal enable_delay   : std_logic             := '0';                   -- Enable signal for delay counter

begin

  -- Pulse generator process: Generates a pulse after a number of clock cycles (based on divisor)
  pulse_generator : process (clock) is
  begin  -- process body
    if rising_edge(clock) then          -- Trigger action on the rising edge of the clock
      if enable_counter = '1' then     -- Only increment the counter if enabled
        counter <= counter + 1;        -- Increment the counter by 1 on each clock cycle
        if counter = divisor then     -- If the counter reaches the divisor
          pulse_out <= '1';           -- Output a pulse to signal baud rate timing
          counter   <= (others => '0'); -- Reset the counter to 0
        else
          pulse_out <= '0';           -- Keep the pulse signal low otherwise
        end if;
      else
        counter <= (others => '0');   -- Reset counter if not enabled
      end if;
    end if;
  end process pulse_generator;

  -- State machine process: Handles state transitions for UART reception (start bit, data bits, stop bit)
  state_machine : process (clock) is
  begin  -- process body
    if rising_edge(clock) then          -- Trigger action on the rising edge of the clock
      case state is
        when idle_s =>  -- Initial idle state: Waiting for start bit (low signal on uart_rx)
          enable_counter <= '0';  -- Disable counter in idle state
          if uart_rx = '0' then  -- When uart_rx is low, start reception process
            state <= start_s;    -- Transition to start state
          end if;
        when start_s =>  -- Start state: Waiting for the start bit to be sampled
          enable_counter <= '1';  -- Enable the counter to start baud rate timing
          if pulse_out = '1' then -- Wait for pulse to indicate baud rate timing
            state <= bit0_s;      -- Transition to receiving bit0
          end if;
        when bit0_s to bit7_s =>  -- States for receiving bits 0 through 7
          if pulse_out = '1' then -- Each time the pulse occurs, sample the next bit
            state <= state'next;   -- Transition to the next bit state
          end if;
        when bit7_s =>  -- Final data bit (bit7)
          if pulse_out = '1' then
            state <= idle_s;      -- After receiving bit7, return to idle state
          end if;
        when others => null;   -- Default case (should not be hit in normal operation)
      end case;
    end if;
  end process state_machine;

  -- Delay line process: Generates a delayed signal for baud rate timing
  delay_line : process (clock) is
  begin  -- process body
    if rising_edge(clock) then          -- Trigger action on the rising edge of the clock
      if pulse_out = '1' then
        -- Start counting once the pulse is output
        enable_delay <= '1';        -- Enable delay counter
      end if;
      
      if delay_counter = half_divisor then  -- If delay counter reaches half the divisor
        enable_delay <= '0';              -- Disable the delay counter
        baudrate_out <= '1';              -- Output the baud rate signal
      else
        baudrate_out <= '0';             -- Keep baud rate signal low otherwise
      end if;
      
      -- If delay is enabled, increment the delay counter
      if enable_delay = '1' then
        delay_counter <= delay_counter + 1;
      else
        delay_counter <= (others => '0');  -- Reset delay counter when not enabled
      end if;
    end if;
  end process delay_line;

end architecture rtl;
