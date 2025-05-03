library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.spike_array.all;

entity top is

  port (

    CLK100MHZ    : in  std_logic;
    uart_txd_in  : in  std_logic;
    uart_rxd_out : out std_logic);

end entity top;

architecture behavior of top is
    -- UART signals
    signal valid_UART      : std_logic;
    signal received_data   : std_logic_vector(7 downto 0);
    -- memory signals
    signal spike_array: spike_array_reg := (others => (others => '0'));
    signal reset_memory    : std_logic;
    signal enable_memory   : std_logic;

    component uart_receiver is
        port (
            clock         : in  std_logic;
            uart_rx       : in  std_logic;
            valid         : out std_logic;
            received_data : out std_logic_vector(7 downto 0)
        );
    end component;

    component memory is
         port (
            valid        : in STD_LOGIC; -- from uart 
            rx           : in  STD_LOGIC_VECTOR(7 downto 0);
            clk          : in  STD_LOGIC;
            reset        : in  STD_LOGIC;
            -- spike_arrays: da 10*N bit, indicizzati downto
            data_reg     : out spike_array_reg; -- inizialmente imposto shift register ad 8 bit per tutte le sinapsi 
            enable       : out STD_LOGIC  
        );
    end component;

begin

        -- Instantiate the UART receiver
    uart_receiver_1 : uart_receiver
        port map (
            clock         => CLK100MHZ,
            uart_rx       => uart_txd_in,
            valid         => valid_UART,
            received_data => received_data
    );

    memory_1 : memory 
         port (
            valid         => valid_UART, 
            rx            => received_data,
            clk           => CLK100MHZ,
            reset         => reset,
            data_reg      => spike_array, -- inizialmente imposto shift register ad 8 bit per tutte le sinapsi 
            enable        => enable_memory 
        );

end architecture behavior;