entity uart_receiver is
  port (
    clock         : in  std_logic;
    uart_rx       : in  std_logic;
    valid         : out std_logic;
    received_data : out std_logic_vector(7 downto 0);
    led           : out std_logic  -- New signal for the LED
  );
end entity uart_receiver;

architecture str of uart_receiver is

  type state_t is (idle_s, start_s, bit0_s, bit1_s, bit2_s, bit3_s, bit4_s, bit5_s, bit6_s, bit7_s, stop_s);
  signal state : state_t := idle_s;

  signal baudrate_out : std_logic;
  signal received_data_s : std_logic_vector(7 downto 0);
  signal sample_wait : boolean := false;

  component sampler_generator is
    port (
      clock        : in  std_logic;
      uart_rx      : in  std_logic;
      baudrate_out : out std_logic
    );
  end component sampler_generator;

begin

  sampler_generator_1 : sampler_generator
    port map (
      clock        => clock,
      uart_rx      => uart_rx,
      baudrate_out => baudrate_out
    );

  main : process (clock) is
  begin
    if rising_edge(clock) then
      case state is
        when idle_s =>
          valid <= '0';
          received_data_s <= (others => '0');
          led <= '0';  -- LED off in idle state
          if uart_rx = '0' then
            -- Detected falling edge (start bit)
            sample_wait <= true;
            state <= start_s;
          end if;

        when start_s =>
          if baudrate_out = '1' then
            if sample_wait = true then
              -- First half-bit delay to sample in middle of start bit
              sample_wait <= false;
            else
              -- Now ready to sample bits
              state <= bit0_s;
            end if;
          end if;

        when bit0_s =>
          if baudrate_out = '1' then
            received_data_s(0) <= uart_rx;
            state <= bit1_s;
          end if;

        when bit1_s =>
          if baudrate_out = '1' then
            received_data_s(1) <= uart_rx;
            state <= bit2_s;
          end if;

        when bit2_s =>
          if baudrate_out = '1' then
            received_data_s(2) <= uart_rx;
            state <= bit3_s;
          end if;

        when bit3_s =>
          if baudrate_out = '1' then
            received_data_s(3) <= uart_rx;
            state <= bit4_s;
          end if;

        when bit4_s =>
          if baudrate_out = '1' then
            received_data_s(4) <= uart_rx;
            state <= bit5_s;
          end if;

        when bit5_s =>
          if baudrate_out = '1' then
            received_data_s(5) <= uart_rx;
            state <= bit6_s;
          end if;

        when bit6_s =>
          if baudrate_out = '1' then
            received_data_s(6) <= uart_rx;
            state <= bit7_s;
          end if;

        when bit7_s =>
          if baudrate_out = '1' then
            received_data_s(7) <= uart_rx;
            state <= stop_s;
          end if;

        when stop_s =>
          if baudrate_out = '1' then
            if uart_rx = '1' then
              -- Good stop bit
              valid <= '1';
              received_data <= received_data_s;
              led <= '1';  -- Turn LED on when the message is valid
            else
              -- Bad stop bit (optional error handling could go here)
              valid <= '0';
              led <= '0';  -- LED off on invalid message
            end if;
            state <= idle_s;
          end if;

        when others =>
          null;
      end case;
    end if;
  end process;

end architecture str;
