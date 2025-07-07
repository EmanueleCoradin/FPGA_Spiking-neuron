library ieee;
use ieee.std_logic_1164.all;  -- Standard logic library, includes std_logic and std_logic_vector
use ieee.numeric_std.all;     -- Numeric library, provides support for unsigned and signed arithmetic

-- Entity declaration: Defines the interface of the baudrate_generator component

entity baudrate_generator is
  port (
    clock        : in  std_logic;         -- Input clock signal
    baudrate_out : out std_logic          -- Output signal for baudrate timing
  );
end entity baudrate_generator;

-- Architecture: Describes the internal implementation of the component (RTL - Register Transfer Level)

architecture rtl of baudrate_generator is
  signal counter   : unsigned(10 downto 0) := (others => '0');       -- 10-bit counter, initialized to 0
  constant divisor : unsigned(10 downto 0) := to_unsigned(867, 11);  -- Constant divisor value for generating the baud rate
begin  -- Start of the RTL architecture

  -- Main process: This process is sensitive to the clock signal (clock)
  
  main : process (clock) is
  begin  -- Process body
    if rising_edge(clock) then            -- Trigger action on the rising edge of the clock signal
      counter <= counter + 1;             -- Increment the counter by 1 on every clock cycle
      
      -- Check if the counter has reached the divisor value
      if counter = divisor then
        baudrate_out <= '1';              -- Output '1' to signal the baud rate timing
        counter      <= (others => '0');  -- Reset the counter to 0 after reaching the divisor value
      else
        baudrate_out <= '0';              -- Keep the baudrate_out signal low until the divisor is reached
      end if;
    end if;
  end process main;                       

end architecture rtl;                    

