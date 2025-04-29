library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use work.spike_array.all;

entity uart_receiver_tb is
end uart_receiver_tb;

architecture behavior of uart_receiver_tb is

  -- Componenti e segnali da testare
  component uart_receiver is
    port (
      valid     : in  std_logic;
      rx        : in  std_logic_vector(7 downto 0);
      clk       : in  std_logic;
      reset     : in  std_logic;
      data_reg  : out spike_array_reg;
      enable    : out std_logic
    );
  end component;

  -- Segnali per test
  signal clk       : std_logic := '0';
  signal reset     : std_logic := '1';
  signal valid     : std_logic := '0';
  signal rx        : std_logic_vector(7 downto 0) := (others => '0');
  signal data_reg  : spike_array_reg;
  signal enable    : std_logic;

  -- Simulazione baudrate generator
  signal baudrate_clk : std_logic := '0';
  constant CLK_PERIOD : time := 10 ns;
  constant BAUD_PERIOD : time := 100 ns;

begin

  -- Clock principale
  clk_process : process
  begin
    while now < 2 ms loop
      clk <= '0';
      wait for CLK_PERIOD / 2;
      clk <= '1';
      wait for CLK_PERIOD / 2;
    end loop;
    wait;
  end process;

  -- Baudrate generator (mock)
  baud_gen_mock : process
  begin
    while now < 2 ms loop
      baudrate_clk <= '0';
      wait for BAUD_PERIOD / 2;
      baudrate_clk <= '1';
      wait for BAUD_PERIOD / 2;
    end loop;
    wait;
  end process;

  -- Collegamento al DUT
  uut: uart_receiver
    port map (
      valid     => valid,
      rx        => rx,
      clk       => baudrate_clk, -- collegato al "baudrate_out" simulato
      reset     => reset,
      data_reg  => data_reg,
      enable    => enable
    );

  -- Stimolo
  stim_proc: process
  begin
    -- Reset iniziale
    wait for 50 ns;
    reset <= '0';

    -- Prima coppia di byte
    wait for BAUD_PERIOD;
    valid <= '1';
    rx <= "10101011"; -- 7 downto 3 = "10101"
    wait for BAUD_PERIOD;

    rx <= "11001100"; -- 7 downto 3 = "11001"
    wait for BAUD_PERIOD;
    valid <= '0';

    -- Seconda coppia di byte
    wait for BAUD_PERIOD * 2;
    valid <= '1';
    rx <= "11110000"; -- 11110
    wait for BAUD_PERIOD;

    rx <= "00011100"; -- 00011
    wait for BAUD_PERIOD;
    valid <= '0';

    wait for BAUD_PERIOD * 5;
    wait;
  end process;

end behavior;

