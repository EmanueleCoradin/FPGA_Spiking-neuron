library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity shift_register is
  generic (N : integer := 16);
  port (
    Clk    : in std_logic;
    Rst    : in std_logic;
    input  : in std_logic;
    output : out std_logic
  );
end shift_register;

architecture Behavioral of shift_register is

  -- Only declare reg if N > 0
  -- Wrap everything in a generate block
begin

  -- Assert right up front to catch bad N values early
  sanity_check : assert N > 0
    report "Generic N must be greater than 0"
    severity failure;

  gen_shift_reg : if N > 0 generate
    signal reg : std_logic_vector(N-1 downto 0);

    begin
      shift_proc : process(Clk, Rst)
      begin
        if Rst = '1' then
          reg <= (others => '0');
        elsif rising_edge(Clk) then
          reg(0) <= input;
          for i in 1 to N-1 loop
            reg(i) <= reg(i - 1);
          end loop;
        end if;
      end process;

      output <= reg(N-1);
  end generate;

end Behavioral;
