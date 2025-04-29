library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity matrix_receiver is
  port (
    clock        : in  std_logic;
    uart_rx      : in  std_logic;
    -- optional outputs for debug:
    debug_valid  : out std_logic;
    debug_row    : out std_logic_vector(9 downto 0)
  );
end entity matrix_receiver;

architecture rtl of matrix_receiver is

  signal received_data : std_logic_vector(7 downto 0);
  signal valid : std_logic;
  
  signal high_byte : std_logic_vector(7 downto 0);
  signal low_byte  : std_logic_vector(7 downto 0);

  type state_t is (wait_first_byte, wait_second_byte);
  signal state : state_t := wait_first_byte;

  -- RAM to store matrix rows
  constant MATRIX_SIZE : integer := 1024; -- Number of rows you want (adjustable)
  type ram_type is array (0 to MATRIX_SIZE-1) of std_logic_vector(9 downto 0);
  signal matrix_mem : ram_type := (others => (others => '0'));

  signal write_address : integer range 0 to MATRIX_SIZE-1 := 0;

  component uart_receiver is
    port (
      clock         : in  std_logic;
      uart_rx       : in  std_logic;
      valid         : out std_logic;
      received_data : out std_logic_vector(7 downto 0)
    );
  end component;

begin

  -- Instantiate UART Receiver
  uart_receiver_inst : uart_receiver
    port map (
      clock         => clock,
      uart_rx       => uart_rx,
      valid         => valid,
      received_data => received_data
    );

  -- Control logic
  process(clock)
    variable temp_row : std_logic_vector(9 downto 0);
  begin
    if rising_edge(clock) then
      debug_valid <= '0';
      case state is
        when wait_first_byte =>
          if valid = '1' then
            high_byte <= received_data;
            state <= wait_second_byte;
          end if;

        when wait_second_byte =>
          if valid = '1' then
            low_byte <= received_data;
            -- Assemble the 10-bit row
            temp_row := high_byte(7 downto 0) & low_byte(7 downto 6); -- 8 bits + 2 bits
            -- Store it in RAM
            matrix_mem(write_address) <= temp_row;
            -- Output for debug
            debug_row <= temp_row;
            debug_valid <= '1';
            -- Increment address
            if write_address < MATRIX_SIZE-1 then
              write_address <= write_address + 1;
            else
              write_address <= 0; -- or stop, depending on your design
            end if;
            -- Return to wait first byte
            state <= wait_first_byte;
          end if;

        when others =>
          state <= wait_first_byte;
      end case;
    end if;
  end process;

end architecture rtl;
