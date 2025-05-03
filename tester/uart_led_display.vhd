library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity uart_led_display is
    port (
        clock           : in  std_logic;                      -- Input clock signal
        uart_rx         : in  std_logic;                      -- UART receive data
        valid           : in  std_logic;                      -- UART valid signal
        received_data   : in  std_logic_vector(7 downto 0);   -- Received data from UART
        led             : out std_logic;                       -- Output to a single LED
        led_enable      : out std_logic                        -- Signal to enable LED output
    );
end entity uart_led_display;

architecture behavior of uart_led_display is

    signal counter       : unsigned(23 downto 0) := (others => '0');   -- 24-bit counter for 1-second delay
    signal bit_index     : integer range 0 to 7 := 0;                  -- Index to track the bit being displayed
    signal message       : std_logic_vector(7 downto 0);                -- Store the received byte
    signal display_ready : boolean := false;                            -- Flag to indicate when to update the LED

    constant clk_freq    : integer := 100000000;  -- Assume 100 MHz clock frequency (adjust as needed)
    constant delay_value : unsigned(23 downto 0) := to_unsigned(clk_freq / 2, 24); -- 1/8 second delay for blinking

begin

    -- Main process to generate the 1/8-second delay and control LED display
    process(clock)
    begin
        if rising_edge(clock) then
            -- Counter to create 1/8-second delay (enough time to blink each bit of the byte)
            if counter = delay_value then
                counter <= (others => '0');  -- Reset counter after 1/8 second
                display_ready <= true;       -- Ready to update the LED
            else
                counter <= counter + 1;     -- Increment counter
                display_ready <= false;      -- Not ready yet
            end if;
        end if;
    end process;

    -- Process to handle the UART received data and display on LED
    process(clock)
    begin
        if rising_edge(clock) then
            if valid = '1' then
                message <= received_data;   -- Store received byte
                bit_index <= 0;             -- Reset bit index when new data arrives
            end if;

            if display_ready then
                -- Display the current bit of the received byte on the LED
                led <= message(bit_index);  -- Display one bit at a time on the LED
                led_enable <= '1';          -- Enable the LED output

                -- Move to the next bit after the current one is displayed
                if bit_index = 7 then
                    bit_index <= 0;        -- Reset to the first bit after displaying the full byte
                else
                    bit_index <= bit_index + 1;
                end if;
            else
                led_enable <= '0';          -- Disable the LED output while not updating
            end if;
        end if;
    end process;

end architecture behavior;
