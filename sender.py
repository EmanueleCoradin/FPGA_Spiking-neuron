import pandas as pd
import numpy as np
import serial
import struct


def encode_for_uart(matrix):
    conversion = np.array([16, 8, 4, 2, 1])
    encoded_data = []

    num_rows, num_columns = matrix.shape
    for col in range(num_columns):
        # Split the column into two parts: first 5 rows and last 5 rows 
        first_part = matrix.iloc[0:5, col].values
        second_part = matrix.iloc[5:10, col].values if num_rows > 5 else []
        # Encode the first part
        first_part_encoded = ''.join(format(int(first_part.dot(conversion)), '05b'))
        
        # Add 3 bits of 0 padding
        if col == num_columns - 1:
            first_part_encoded += '111'
        else:
            first_part_encoded += '000'
        # Encode the second part
        second_part_encoded = ''.join(format(int(second_part.dot(conversion)), '05b'))
        # Add 3 bits of 0 padding
        if col == num_columns - 1:
            second_part_encoded += '111'
        else:
            second_part_encoded += '000'

        # Concatenate both parts
        full_encoded = first_part_encoded + second_part_encoded

        # Store the encoded data for this column
        encoded_data.append(full_encoded)

    return np.array(encoded_data)

# Open UART port
ser = serial.Serial('/dev/ttyUSB0', baudrate=115200)
print(f"Connection opened.")

# Load the matrix from the CSV file
file_path = './CMS-SpikingNeuralNetwork/matrix.csv'
matrix = pd.read_csv(file_path, header=None, nrows=10)

print(f"Encoding data.")
# Encode the matrix
encoded_matrix = encode_for_uart(matrix)
print(f"Sending data.")
# Display a small portion of the encoded data
for col in encoded_matrix:
    # Extract the column (as a byte array)
    col_data = col.tobytes()
    
    # Send the column data via UART
    ser.write(col_data)
print(f"Data sent via UART.")
ser.close()
