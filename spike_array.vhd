library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

-- Ricevo in input dalla UART un array a 8 bit: i primi 4 rappresentano l'indirizzo
-- della sinapsi, mentre il successivo è il bit dell'informazione da trasmettere 
package spike_array is 
    type spike_array_reg is array (0 to 9) of STD_LOGIC_VECTOR(1024-1 to 0);
end package;
