import pandas as pd
import numpy as np
import serial
import struct
import time


def encode_for_uart(matrix):
	conversion = np.array([16, 8, 4, 2, 1])
	encoded_data = []

	num_rows, num_columns = matrix.shape
	for col in range(num_columns):
		# Split the column into two parts: first 5 rows and last 5 rows 
		first_part = matrix.iloc[0:5, col].values
		second_part = matrix.iloc[5:10, col].values if num_rows > 5 else []
		# Encode the first part
		first_part_encoded = int(first_part.dot(conversion))

		# Add 3 bits of 0 padding
		if col == num_columns - 1:
			first_part_encoded = first_part_encoded*8+7
		else:
			first_part_encoded += first_part_encoded*8
		first_part_encoded = first_part_encoded.to_bytes(1, 'big')
		encoded_data.append(first_part_encoded)

		# Encode the second part
		second_part_encoded = int(second_part.dot(conversion))
		# Add 3 bits of 0 padding
		if col == num_columns - 1:
			second_part_encoded = second_part_encoded*8+7
		else:
			second_part_encoded += second_part_encoded*8
		second_part_encoded = second_part_encoded.to_bytes(1, 'big')

		encoded_data.append(second_part_encoded)

	return encoded_data

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
for col in [encoded_matrix[0]]:
	print(col)

	# Send the column data via UART
	ser.write(col)

print(f"Data sent via UART.")
ser.close()