library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity output_memory2 is
  generic (
    DEPTH        : integer;          -- Number of memory columns (samples)
    DELAY_CYCLES : integer      -- Delay cycles before transmission
  );
  port (
    rst           : in  std_logic;                                  -- Reset
    clk           : in  std_logic;                                  -- Clock
    UART_busy     : in  std_logic;                                  -- UART busy signal
    incoming_data : in  std_logic_vector(7 downto 0);               -- 8-bit voltage input
    valid_UART    : out std_logic;                                  -- UART transmit trigger (1-cycle pulse)
    outgoing_data : out std_logic_vector(7 downto 0);               -- Output voltage for UART
    write_trigger : in  std_logic                                   -- Write enable trigger
  );
end output_memory2;



architecture Behavioral of output_memory2 is

  -- Memory: 8 rows (channels), each with DEPTH samples
  type voltage_reg is array (7 downto 0) of std_logic_vector(DEPTH-1 downto 0);
      -- Register for all the voltage signals. Collects DEPTH points in an 8 x DEPTH matrix, where each column represents an instance of signal
    -- and each row represents a time instant
    
    -- I will use the terms "rows" and "columns" quite liberally and sometimes when I visualize stuff in my head, people tell me that 
    -- I don't make a lot of sense.
    -- As a reference, this is how I envision the memory:
    

    --  0 1 2 3 4 5 6 7 8 9  ... DEPTH-2 DEPTH-1
    --  _ _ _ _ _ _ _ _ _ _  ...  _ _
    -- | | | | | | | | | | | ... | | |   row 0
    --  - - - - - - - - - -  ...  - -  
    -- | | | | | | | | | | | ... | | |   row 1
    --  - - - - - - - - - -  ...  - -
    -- ... ... ... ... ... ... ... ...
    --  _ _ _ _ _ _ _ _ _ _  ...  _ _
    -- | | | | | | | | | | | ... | | |   row 7
    --  - - - - - - - - - -  ...  - -
    
  signal voltage_register : voltage_reg := (others => (others => '1'));

  signal column_index      : integer range 0 to DEPTH-1 := DEPTH-1;
  signal read_index        : integer range 0 to DEPTH-1 := DEPTH-1;
  signal memory_full       : std_logic := '0';
  signal currently_sending : std_logic := '0';
  signal UART_pulse        : std_logic := '0';
  signal outgoing_buffer   : std_logic_vector(7 downto 0);
  signal delay_counter     : integer range 0 to DELAY_CYCLES := 0;
  signal delay_done        : std_logic := '0';
  signal uart_busy_prev : std_logic := '0';

begin

  -- Output assignments
  valid_UART    <= UART_pulse;
  outgoing_data <= outgoing_buffer;

  process(clk, rst)
  begin
    if rst = '1' then
      voltage_register   <= (others => (others => '0'));
      column_index       <= DEPTH - 1;
      read_index         <= DEPTH - 1;
      memory_full        <= '0';
      currently_sending  <= '0';
      UART_pulse         <= '0';
      delay_counter      <= 0;
      delay_done         <= '0';

    elsif rising_edge(clk) then
      UART_pulse <= '0';

      -- Phase 1: Write incoming data
      if memory_full = '0' then
        if write_trigger = '1' then
          for i in 0 to 7 loop
            voltage_register(i)(column_index) <= incoming_data(i);
          end loop;

          if column_index = 0 then
            column_index       <= DEPTH - 1;
            read_index         <= DEPTH - 1;
            memory_full        <= '1';
            delay_counter      <= 0;
            delay_done         <= '0';
            currently_sending  <= '0';
          else
            column_index <= column_index - 1;
          end if;
        end if;
      end if;

      -- Phase 2: Optional delay before sending
      if memory_full = '1' and delay_done = '0' then
        if delay_counter < DELAY_CYCLES - 1 then
          delay_counter <= delay_counter + 1;
        else
          delay_done        <= '1';
          currently_sending <= '1';
        end if;
      end if;
    
     -- Phase 3: Send data via UART
    if currently_sending = '1' then
      -- Detect falling edge of UART_busy
      if (uart_busy_prev = '1' and UART_busy = '0') or read_index = DEPTH - 1 then
        for i in 0 to 7 loop
          outgoing_buffer(i) <= voltage_register(i)(read_index);
        end loop;
    
        UART_pulse <= '1'; -- Pulse valid_UART for 1 clock cycle
    
        if read_index = 0 then
          currently_sending <= '0'; -- Transmission complete
        else
          read_index <= read_index - 1;
        end if;
      end if;
    end if;
    
    -- Track previous UART_busy state
    uart_busy_prev <= UART_busy;
  end if;
  end process;

end Behavioral;
