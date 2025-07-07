library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity uart_receiver is
  port (
    clk           : in  std_logic;                            -- 100 MHz clock
    uart_rx       : in  std_logic;                            -- UART RX input
    received_data : out std_logic_vector(7 downto 0);         -- 8-bit output data
    valid         : out std_logic                             -- Pulse when byte is received
    --error       : out std_logic                             -- Optional: Stop bit error (not used)
  );
end entity;

architecture rtl of uart_receiver is

  -- FSM states
  type state_t is (idle, start, data, stop);
  signal state : state_t := idle;

  -- UART timing parameters
  constant oversample   : integer := 16;                      -- 16x oversampling
  constant baud_div     : integer := 867;                     -- 100 MHz / 115200 ≈ 867
  constant sample_div   : integer := baud_div / oversample;   -- ≈ 54

  -- Timing counters
  signal sample_cnt     : integer range 0 to sample_div-1 := 0;
  signal oversample_cnt : integer range 0 to oversample-1 := 0;

  -- Data reception
  signal bit_index      : integer range 0 to 7 := 0;
  signal rx_shift       : std_logic_vector(7 downto 0) := (others => '0');

  -- RX input synchronization
  signal rx_sync        : std_logic_vector(1 downto 0) := (others => '1');
  signal rx_filtered    : std_logic;

  -- Pulse indicating when to sample
  signal sample_tick    : std_logic := '0';

begin

  -------------------------------------------------------------------
  -- Synchronize asynchronous RX input to the system clock domain
  -------------------------------------------------------------------
  process(clk)
  begin
    if rising_edge(clk) then
      rx_sync <= rx_sync(0) & uart_rx;
    end if;
  end process;

  rx_filtered <= rx_sync(1);  -- Debounced and synchronized RX line

  -------------------------------------------------------------------
  -- Generate a sampling tick every sample_div clock cycles
  -------------------------------------------------------------------
  process(clk)
  begin
    if rising_edge(clk) then
      if sample_cnt = sample_div - 1 then
        sample_cnt  <= 0;
        sample_tick <= '1';
      else
        sample_cnt  <= sample_cnt + 1;
        sample_tick <= '0';
      end if;
    end if;
  end process;

  -------------------------------------------------------------------
  -- Main FSM: Receives UART frame (1 start + 8 data + 1 stop)
  -------------------------------------------------------------------
  process(clk)
  begin
    if rising_edge(clk) then
      valid <= '0';  -- Default: No valid data this cycle
      --error <= '0'; -- Uncomment if stop bit error is used

      if sample_tick = '1' then
        case state is

          -----------------------------------------------------------
          -- Wait for start bit (RX goes low)
          -----------------------------------------------------------
          when idle =>
            if rx_filtered = '0' then
              oversample_cnt <= 0;
              state <= start;
            end if;

          -----------------------------------------------------------
          -- Sample start bit at midpoint (oversample count = 7)
          -----------------------------------------------------------
          when start =>
            if oversample_cnt = 7 then
              if rx_filtered = '0' then  -- Valid start bit
                oversample_cnt <= 0;
                bit_index <= 0;
                state <= data;
              else                       -- False start bit
                state <= idle;
              end if;
            else
              oversample_cnt <= oversample_cnt + 1;
            end if;

          -----------------------------------------------------------
          -- Receive 8 data bits, sampled at center of each bit period
          -----------------------------------------------------------
          when data =>
            if oversample_cnt = 15 then
              rx_shift(bit_index) <= rx_filtered;
              oversample_cnt <= 0;
              if bit_index = 7 then
                state <= stop;
              else
                bit_index <= bit_index + 1;
              end if;
            else
              oversample_cnt <= oversample_cnt + 1;
            end if;

          -----------------------------------------------------------
          -- Sample stop bit, check for validity
          -----------------------------------------------------------
          when stop =>
            if oversample_cnt = 15 then
              if rx_filtered = '1' then       -- Valid stop bit
                received_data <= rx_shift;
                valid <= '1';
              else                            -- Framing error (optional)
                valid <= '1';
                --error <= '1';
              end if;
              state <= idle;
              oversample_cnt <= 0;
            else
              oversample_cnt <= oversample_cnt + 1;
            end if;

        end case;
      end if;
    end if;
  end process;

end architecture;
