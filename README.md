# FPGA Spiking Neuron

This project implements a hardware spiking neuron model using VHDL, designed for FPGA. It includes tools for generating and sending stimulus data to the design.

## 📦 Project Structure

FPGA_Spiking-neuron/
├── stimuli_files/ # Input stimulus data
├── stimuli_generator/ # Code to generate stimuli (C++)
├── vhdl_sources/ # VHDL source files (spiking neuron logic)
├── sender.py # Python script to send stimuli (e.g., via serial or socket)
├── README.md # You're here


## 🔬 Overview

This repository contains:

- A **VHDL implementation** of a spiking neuron model
- A **C++ stimulus generator**
- A **Python sender script** to transmit input data to the FPGA board via UART
- Example stimulus files for testing and simulation

The VHDL model has been tested on a Xilinx Artix-7 FPGA board.

## 🧠 Neuron Model

The current design implements a simplified Leaky Integrate-and-Fire (LIF)] spiking neuron model. The neuron integrates input stimuli over time and emits a spike when the membrane potential exceeds a threshold.

## 🚀 Getting Started

### Prerequisites

- **FPGA toolchain**: e.g., Xilinx Vivado
- **C++ compiler**: For the stimulus generator
- **Python 3.x**: For sending stimuli. 

Python scripts require the following packages:

- `numpy`
- `pandas`
- `matplotlib`
- `pyserial`

CERN ROOT framework is required as well. 
