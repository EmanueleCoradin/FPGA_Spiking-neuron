library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity top_uart_display is
    port (
        clock           : in  std_logic;                      -- Input clock signal
        uart_rx         : in  std_logic;                      -- UART receive data
        led             : out  std_logic_vector(7 downto 0)   -- Output to a single LED
    );
end entity top_uart_display;

architecture behavior of top_uart_display is

    signal valid           : std_logic;
    signal received_data   : std_logic_vector(7 downto 0);

    -- ILA (Integrated Logic Analyzer) probe signals
    signal ila_received_data : std_logic_vector(7 downto 0);  -- Probe signal for received data

    -- UART Receiver component declaration
    component uart_receiver is
        port (
            clock         : in  std_logic;
            uart_rx       : in  std_logic;
            valid         : out std_logic;
            received_data : out std_logic_vector(7 downto 0)
        );
    end component;

    -- ILA component declaration
    --component ila is
    --    port (
    --        clk      : in  std_logic;                 -- Clock signal
    --        probe0   : in  std_logic_vector(7 downto 0);  -- Probe signal for received_data
    --        probe1   : in  std_logic  -- Probe signal for received_data
    --    );
    --end component;

begin

    -- Instantiate the UART receiver
    uart_receiver_1 : uart_receiver
        port map (
            clock         => clock,
            uart_rx       => uart_rx,
            valid         => valid,
            received_data => received_data
        );

    -- Instantiate the ILA to probe the received_data signal
    --ila_1: ila
    --    port map (
    --        clk      => clock,                    -- Provide the clock to the ILA
    --        probe0   => received_data, -- Data signal you want to monitor
    --       probe1 => valid          -- Valid signal to trigger when data is valid
    --    );

    -- Process to update LED based on valid data
    process(clock)
    begin
        if rising_edge(clock) then
            if valid = '1' then
                led <= received_data;  -- Output received data to LEDs when valid
            end if;
        end if;
    end process;

end architecture behavior;

