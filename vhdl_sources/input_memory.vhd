library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity memory is
  generic (
    N            : integer;  -- Number of memory slots
    INDEX_WIDTH  : integer   -- Width of read_index (must satisfy 2^INDEX_WIDTH ≥ N)
  );
  port (
    valid      : in  STD_LOGIC;
    rx         : in  STD_LOGIC_VECTOR(7 downto 0);  -- Incoming byte from UART
    clk        : in  STD_LOGIC;
    reset      : in  STD_LOGIC;
    out_spike  : out STD_LOGIC_VECTOR(9 downto 0);  -- One-hot spikes 
    enable     : out STD_LOGIC;                     -- Indicates memory is ready
    button     : in  STD_LOGIC                      -- Trigger memory read mode
  );
end entity memory;

architecture Behavioral of memory is

  -- Memory for 10 input spike streams, each N bits deep
  type spike_array_reg is array (0 to 9) of STD_LOGIC_VECTOR(N-1 downto 0);
  signal out_spike_array : spike_array_reg := (others => (others => '0'));

  -- Control signals
  signal waiting_second_byte_internal : std_logic := '0';
  signal write_position_int           : integer range 0 to N-1 := N-1;
  signal memory_ready                 : std_logic := '0';
  signal set_memory_ready             : std_logic := '0';
  signal read_index_int               : integer range 0 to N-1 := N-1;
  signal write_done                   : std_logic := '0';

  -- Output buffer for one spike vector
  signal out_spike_bit : std_logic_vector(9 downto 0);

  -- Button debouncing
  signal button_debounced : std_logic := '0';
  signal button_sync_0    : std_logic := '0';
  signal button_sync_1    : std_logic := '0';
  signal button_cnt       : integer range 0 to 1000000 := 0;
  constant DEBOUNCE_LIMIT : integer := 1000000;  -- ~10ms @100MHz

  -- Edge detection for valid signal
  signal valid_prev : std_logic := '0';
  signal valid_edge : std_logic := '0';

  -- Synthesis hint to keep signals from being optimized away
  attribute keep : boolean;
  attribute keep of write_position_int : signal is true;
  attribute keep of out_spike_array    : signal is true;

begin

  -- Debouncer for button input
  process(clk, reset)
  begin
    if reset = '1' then
      button_sync_0     <= '0';
      button_sync_1     <= '0';
      button_cnt        <= 0;
      button_debounced  <= '0';
    elsif rising_edge(clk) then
      button_sync_0 <= button;
      button_sync_1 <= button_sync_0;

      if button_sync_1 = '1' then
        if button_cnt < DEBOUNCE_LIMIT then
          button_cnt <= button_cnt + 1;
        end if;
      else
        button_cnt <= 0;
      end if;

      if button_cnt = DEBOUNCE_LIMIT then
        button_debounced <= '1';
      else
        button_debounced <= '0';
      end if;
    end if;
  end process;

  -- Main memory logic
  process(clk, reset)
  begin
    if reset = '1' then
      out_spike_array              <= (others => (others => '0'));
      waiting_second_byte_internal <= '0';
      write_position_int           <= N-1;
      memory_ready                 <= '0';
      set_memory_ready             <= '0';
      write_done                   <= '0';
      out_spike_bit                <= (others => '0');
      valid_prev                   <= '0';

    elsif rising_edge(clk) then

      -- Detect rising edge of valid signal
      if valid = '1' and valid_prev = '0' then
        valid_edge <= '1';
      else
        valid_edge <= '0';
      end if;
      valid_prev <= valid;

      -- Output spikes sequentially after memory is filled
      if write_done = '1' and read_index_int > 0 then
        for i in 0 to 9 loop
          out_spike_bit(i) <= out_spike_array(i)(read_index_int);
        end loop;
        read_index_int <= read_index_int - 1;
      else
        out_spike_bit <= (others => '0');
      end if;

      -- Writing data into memory (two-byte protocol)
      if valid_edge = '1' and memory_ready = '0' and write_done = '0' then
        if waiting_second_byte_internal = '0' then
          -- First byte → spike bits 0-4 from rx(3) to rx(7)
          for i in 0 to 4 loop
            out_spike_array(i)(write_position_int) <= rx(i + 3);
          end loop;
          waiting_second_byte_internal <= '1';

        else
          -- Second byte → spike bits 5-9 from rx(3) to rx(7)
          for i in 5 to 9 loop
            out_spike_array(i)(write_position_int) <= rx(i - 2);
          end loop;
          waiting_second_byte_internal <= '0';
          
          -- Decrement write pointer or wrap around
          if write_position_int = 0 then
            write_position_int <= N-1;
          else
            write_position_int <= write_position_int - 1;
          end if;
          -- If second byte indicates end of transmission, prepare to output
          if rx(2 downto 0) = "111" then
            set_memory_ready <= '1';
          end if;
        end if;

      elsif memory_ready = '1' then
        memory_ready <= '0';  -- Clear ready flag after one cycle
      end if;

      -- Button press triggers transition to output mode
      if button_debounced = '1' and  write_done='0' then
        set_memory_ready <= '1';
      end if;

      -- Set memory ready and begin output sequence
      if set_memory_ready = '1' then
        set_memory_ready   <= '0';
        memory_ready       <= '1';
        write_done         <= '1';
        write_position_int <= N-1;
      end if;

    end if;
  end process;

  -- Assign outputs
  out_spike <= out_spike_bit;
  enable    <= memory_ready;

end Behavioral;
