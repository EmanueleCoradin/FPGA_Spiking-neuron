library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity top is
  generic (
    top_N            : integer := 1024;      -- Number of memory slots
    top_INDEX_WIDTH  : integer := 9;         -- Width of read_index (must satisfy 2^INDEX_WIDTH ≥ N)
    top_DEPTH        : integer := 1024;      -- Number of memory columns (samples)
    top_DELAY_CYCLES : integer := 1_000_000  -- Delay cycles before transmission
  );
  port (
    CLK100MHZ    : in  std_logic;
    uart_txd_in  : in  std_logic;
    button       : in  std_logic;
    led          : out std_logic_vector(2 downto 0);  
    uart_tx      : out std_logic                       -- UART TX output
  );
end entity top;

architecture behavior of top is

  -- UART reception signals
  signal valid_UART        : std_logic;
  signal received_data_s   : std_logic_vector(7 downto 0);
  signal button_s          : std_logic;

  -- Memory control signals
  signal reset_memory      : std_logic := '0';
  signal enable_memory     : std_logic := '1';
  signal out_spike         : std_logic_vector(9 downto 0);  -- Output vector of spikes from memory
  signal enable            : std_logic;                     -- Declared but unused
  signal start_simulation  : boolean := FALSE;              -- Boolean used to track simulation start

  -- SNN input and output signals
  signal reset_n           : std_logic := '0';
  signal sp_0              : std_logic := '0';
  signal sp_1              : std_logic := '0';
  signal sp_2              : std_logic := '0';
  signal sp_3              : std_logic := '0';
  signal sp_4              : std_logic := '0';
  signal sp_5              : std_logic := '0';
  signal sp_6              : std_logic := '0';
  signal sp_7              : std_logic := '0';
  signal sp_8              : std_logic := '0';
  signal sp_9              : std_logic := '0';
  signal spike_out_0       : std_logic;                    -- Declared but not connected to anything
  signal voltage_out_0     : std_logic_vector(7 downto 0);

  -- Interconnect signals between memory, SNN, and UART
  signal rst               : std_logic := '0';             -- Reset input for output memory
  signal UART_busy         : std_logic;
  signal valid_UART_trx    : std_logic;
  signal outgoing_data     : std_logic_vector(7 downto 0);

  -- Write trigger signal used to initiate write in output memory
  signal write_trigger     : std_logic := '0';

  -- UART receiver component declaration
  component uart_receiver is
    port (
      clk           : in  std_logic;
      uart_rx       : in  std_logic;
      valid         : out std_logic;
      received_data : out std_logic_vector(7 downto 0)
    );
  end component;

  -- Memory component receives UART data and outputs spike vector
  component memory is
    generic (
      N           : integer;
      INDEX_WIDTH : integer
    );
    port (
      valid      : in  std_logic;
      rx         : in  std_logic_vector(7 downto 0);
      clk        : in  std_logic;
      reset      : in  std_logic;
      out_spike  : out std_logic_vector(9 downto 0);
      enable     : out std_logic;
      button     : in  std_logic
    );
  end component;

  -- SNN core component
  component snn_architecture is
    port (
      clk           : in  std_logic;
      reset_n       : in  std_logic;
      sp_0          : in  std_logic;
      sp_1          : in  std_logic;
      sp_2          : in  std_logic;
      sp_3          : in  std_logic;
      sp_4          : in  std_logic;
      sp_5          : in  std_logic;
      sp_6          : in  std_logic;
      sp_7          : in  std_logic;
      sp_8          : in  std_logic;
      sp_9          : in  std_logic;
      spike_out_0   : out std_logic;
      voltage_out_0 : out std_logic_vector(7 downto 0)
    );
  end component;

  -- Output memory stores and serializes SNN output over time
  component output_memory2 is
    generic (
      DEPTH        : integer;
      DELAY_CYCLES : integer
    );
    port (
      rst           : in  std_logic;
      clk           : in  std_logic;
      UART_busy     : in  std_logic;
      incoming_data : in  std_logic_vector(7 downto 0);
      valid_UART    : out std_logic;
      outgoing_data : out std_logic_vector(7 downto 0);
      write_trigger : in  std_logic
    );
  end component;

  -- UART transmitter sends bytes out over serial line
  component uart_transmitter is
    port (
      clock        : in  std_logic;
      data_to_send : in  std_logic_vector(7 downto 0);
      data_valid   : in  std_logic;
      busy         : out std_logic;
      uart_tx      : out std_logic
    );
  end component;

begin

  -- UART receiver instantiation
  uart_receiver_1 : uart_receiver
    port map (
      clk           => CLK100MHZ,
      uart_rx       => uart_txd_in,
      valid         => valid_UART,
      received_data => received_data_s
    );

  -- Memory instantiation: receives UART data and outputs spike pattern
  memory_1 : memory
    generic map (
      N           => top_N,
      INDEX_WIDTH => top_INDEX_WIDTH
    )
    port map (
      valid      => valid_UART,
      rx         => received_data_s,
      clk        => CLK100MHZ,
      reset      => reset_memory,
      out_spike  => out_spike,
      enable     => enable_memory,
      button     => button_s
    );

  -- SNN module receives spike inputs and produces voltage/spike outputs
  snn_architecture_1 : snn_architecture
    port map (
      clk           => CLK100MHZ,
      reset_n       => reset_n,
      sp_0          => sp_0,
      sp_1          => sp_1,
      sp_2          => sp_2,
      sp_3          => sp_3,
      sp_4          => sp_4,
      sp_5          => sp_5,
      sp_6          => sp_6,
      sp_7          => sp_7,
      sp_8          => sp_8,
      sp_9          => sp_9,
      spike_out_0   => spike_out_0,           -- Output declared, but not used further
      voltage_out_0 => voltage_out_0
    );

  -- Output memory buffers voltage outputs and schedules UART transmission
  output_mem_inst : output_memory2
    generic map (
      DEPTH        => top_DEPTH,
      DELAY_CYCLES => top_DELAY_CYCLES
    )
    port map (
      rst           => rst,
      clk           => CLK100MHZ,
      UART_busy     => UART_busy,
      incoming_data => voltage_out_0,
      valid_UART    => valid_UART_trx,
      outgoing_data => outgoing_data,
      write_trigger => write_trigger
    );

  -- UART transmitter sends buffered output from memory
  uart_tx_inst : uart_transmitter
    port map (
      clock        => CLK100MHZ,
      data_valid   => valid_UART_trx,
      data_to_send => outgoing_data,
      busy         => UART_busy,
      uart_tx      => uart_tx
    );

  -- Process that controls the simulation sequence:
  -- 1. Wait for memory enable
  -- 2. Trigger spike transfer to SNN
  -- 3. Trigger output memory write
  process(CLK100MHZ)
  begin
    if rising_edge(CLK100MHZ) then
      if enable_memory = '1' then
        led(1)           <= '1';
        reset_n          <= '0';    -- Assert reset on SNN
        start_simulation <= TRUE;
        write_trigger    <= '1';    -- Start writing to output memory

      elsif not start_simulation then
        led(0)          <= '1';
        led(1)          <= '0';
        led(2)          <= '0';
        write_trigger   <= '0';    -- Idle state
        reset_n         <= '1';    -- Hold SNN in reset
        -- Clear all spike inputs
        sp_0 <= '0'; sp_1 <= '0'; sp_2 <= '0'; sp_3 <= '0'; sp_4 <= '0';
        sp_5 <= '0'; sp_6 <= '0'; sp_7 <= '0'; sp_8 <= '0'; sp_9 <= '0';

      elsif start_simulation then
        led(0) <= '0';
        led(2) <= '1';
        -- Assign spikes from memory output to SNN
        sp_0 <= out_spike(0);
        sp_1 <= out_spike(1);
        sp_2 <= out_spike(2);
        sp_3 <= out_spike(3);
        sp_4 <= out_spike(4);
        sp_5 <= out_spike(5);
        sp_6 <= out_spike(6);
        sp_7 <= out_spike(7);
        sp_8 <= out_spike(8);
        sp_9 <= out_spike(9);
      end if;
    end if;
  end process;

  -- Button passthrough (to debounced version or internal logic)
  button_s <= button;

end architecture behavior;
