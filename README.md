# Identifying Network Protocols using Neural Networks for IoT Applications
I experiment with various neural network architectures to build a model that can effectively and efficiently identify network protocols (BLE and Zigbee) and distinguish them from noise given a raw RF signal (in I/Q data form). I find that FCNs (fully convolutional networks) have the best tradeoff between accuracy, performance, and inference on arbitrary number of samples.

My report PDF can be found [here](https://github.com/eminalparslan/protocol-identification/blob/main/project_writeup_emin_arslan.pdf).

The model I ended up selected based on the results of the experiments was the regularized FCN (fully convolutional network) model. My goal was to make the model big enough to the point where it could be useful in practice, but not too big such that it would be difficult to run on IoT devices.
