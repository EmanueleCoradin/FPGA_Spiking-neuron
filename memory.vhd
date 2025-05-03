library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use work.spike_array.all;

entity memory is
  -- generic( N : integer := 1000 );
  port (
    valid : in STD_LOGIC; -- dalla uart 
    rx           : in  STD_LOGIC_VECTOR(7 downto 0);
    clk          : in  STD_LOGIC;
    reset        : in  STD_LOGIC;
    -- spike_arrays: da 10*N bit, indicizzati downto
    data_reg  : out spike_array_reg; -- inizialmente imposto shift register ad 8 bit per tutte le sinapsi 
    enable : out STD_LOGIC  
  ); -- mando in output gli array con i LSB finali per ogni sinapsi 
end entity memory;

architecture Behavioral of memory is

  signal out_spike_array: spike_array_reg := (others => (others => '0')); -- all'inizio metto a zero tutti gli elementi 
  signal waiting_second_byte : boolean := FALSE; 
  signal first_half : STD_LOGIC_VECTOR(7 downto 0);
  signal write_position : integer range 0 to 1023 := 1023;
  signal spike_data : STD_LOGIC_VECTOR(9 downto 0);
  signal memory_ready:  STD_LOGIC:='0';


begin
  process(clk, reset)
  begin
    if reset = '1' then
       out_spike_array <= (others => (others => '0'));
       waiting_second_byte <= FALSE;
       write_position <= 1023;
       memory_ready <= '0';
       
    elsif rising_edge(clk) then
      -- Esegue logica con la frequenza di clk
      if valid = '1' and memory_ready = '0' then 
        if not waiting_second_byte then 
          -- prima ricezione della stringa 
          first_half <= rx;
          waiting_second_byte <= TRUE;
        else 
          -- seconda ricezione 
          spike_data <= first_half(7 downto 3) & rx(7 downto 3);

          -- scrivo nella memoria
          for i in 0 to 9 loop
            out_spike_array(i)(write_position) <= spike_data(i);
          end loop;
 
          -- protezione da underflow 
          if write_position = 0 then 
            write_position <= 1023;
          else
            write_position <= write_position - 1;
          end if; 

          waiting_second_byte <= FALSE;
          
          if rx(2 downto 0) = "111" then
            memory_ready <= '1';
          end if;
        end if;
      elsif memory_ready = '1' then
        memory_ready <= '0';
      end if;

    end if;
  end process;
 
data_reg <= out_spike_array; -- invio in output lo spike array con l'enable 
enable <= memory_ready; -- quando attivare l'enable 
                        
end Behavioral;


