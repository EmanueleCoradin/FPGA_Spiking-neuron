import pandas as pd
import numpy as np
import serial
import struct
import time
import matplotlib.pyplot as plt
import threading
import sys

# === Global flag to signal when to stop UART reading ===
stop_reading = False

def wait_for_stop():
    """
    Waits for the user to press ENTER to set the stop_reading flag,
    signaling the main thread to stop reading from UART.
    """
    global stop_reading
    input("Press ENTER to stop reading from UART...\n")
    stop_reading = True

def encode_for_uart(matrix):
    """
    Encodes a 10xN pandas DataFrame of binary values (0 or 1) into UART-friendly bytes.

    The matrix is split into two parts per column (rows 0-4 and 5-9),
    each part is reversed, dot-multiplied with a fixed weight vector,
    then scaled and adjusted for the last column.

    Args:
        matrix (pd.DataFrame): Input binary matrix (10 rows, N columns).

    Returns:
        List of bytes: Encoded data bytes ready to be sent over UART.
    """
    conversion = np.array([16, 8, 4, 2, 1])  # weights for binary encoding (bit significance)
    encoded_data = []

    num_rows, num_columns = matrix.shape

    for col in range(num_columns):
        # Reverse the order of rows in each 5-row group for correct bit mapping
        first_part = matrix.iloc[0:5, col].values[::-1]  
        second_part = matrix.iloc[5:10, col].values[::-1] if num_rows > 5 else [] 

        try:
            # Dot product maps the 5-bit vector into an integer
            first_part_encoded = int(first_part.dot(conversion))
            # Multiply by 8 (shift left by 3 bits), and add 7 if last column for protocol marker
            first_part_encoded = first_part_encoded * 8 + (7 if col == num_columns - 1 else 0)
        except Exception as e:
            print(f"Error encoding first part at column {col}: {first_part}, Exception: {e}")
            sys.exit(1)

        try:
            # Same encoding for second part, or zero if empty
            second_part_encoded = int(second_part.dot(conversion)) if len(second_part) > 0 else 0
            second_part_encoded = second_part_encoded * 8 + (7 if col == num_columns - 1 else 0)
        except Exception as e:
            print(f"Error encoding second part at column {col}: {second_part}, Exception: {e}")
            sys.exit(1)

        # Append encoded bytes to list (big endian 1 byte each)
        encoded_data.append(first_part_encoded.to_bytes(1, 'big'))
        encoded_data.append(second_part_encoded.to_bytes(1, 'big'))

    return encoded_data

# === UART SETUP ===
ser = serial.Serial('/dev/ttyUSB1', baudrate=115200, timeout=10)
print("UART connection opened.")

# === LOAD AND ENCODE MATRIX DATA ===
matrix = pd.read_csv('matrix.csv', header=None, nrows=10)
print("Encoding data...")
encoded_matrix = encode_for_uart(matrix)

# === SEND ENCODED DATA TO FPGA VIA UART ===
print("Sending data...")
for index, col_bytes in enumerate(encoded_matrix):
    print(f"Sending byte {index}: {col_bytes}")
    ser.write(col_bytes)

# === START THREAD TO WAIT FOR USER TO STOP READING ===
threading.Thread(target=wait_for_stop, daemon=True).start()

# === CONTINUOUSLY READ DATA FROM UART UNTIL USER STOPS ===
print("Reading data... Press ENTER to stop.")
received_data = bytearray()

while not stop_reading:
    if ser.in_waiting:
        # Read all available bytes at once and append to buffer
        received_data.extend(ser.read(ser.in_waiting))

print(f"\nStopped reading. Total bytes received: {len(received_data)}")
print(received_data)

# === CONVERT RECEIVED BYTE DATA TO SIGNED INTEGERS ===
int_data = np.frombuffer(received_data, dtype=np.int8)  # interpret bytes as signed 8-bit integers
df = pd.DataFrame(int_data, columns=["value"])

# === PLOT THE RECEIVED DATA ===
plt.figure(figsize=(12, 5))
plt.plot(df["value"], label="UART Output (signed)")
plt.title("UART Received Data (signed interpretation)")
plt.xlabel("Byte Index")
plt.ylabel("Value (-128 to 127)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# === CLOSE UART CONNECTION ===
ser.close()
print("UART connection closed.")

