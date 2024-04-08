# Eye Blink Detection on EEG Data
An algorithm to detect eye blinks on a low channel EEG recording without any supervised learning. 

#### Project Structure

*  **detection_algorithm:** contains the algorithm for the detection
*  **sample_data:** contains a sample eeg recording

#### Running the project

```shell
python -m venv env
```
```shell
source env/bin/activate
```

```shell
pip install -r requirements.txt
```

```shell
cd detection_algorithm
```
```shell
python eyeblink_detection.py
```



#### Raw EEG data
Data is stored in a .csv format where column 0 represents time, and column 1 and 2 represents raw EEG potentials (in uV) for channel 1 (Fp1) and channel 2 (Fp2) respectively.

#### Output
* The start and end time of each eyeblink is returned by the algorithm 

#### Constraints
The goal of the project is to introduce BCI into consumer hardware. In order to serve that purpose, the project was created with the following constraints in mind
* Low channel EEG recording. The algorithm is generalizable to multiple channels but can run on as low as a 2 channel recording. 
* Unsupervised training. We don't need to retrain and calibrate the BCI hardware to each user. The algorithm is generalizable to any user without any user calibration.


#### Usecases
* The algorithm can be used with EEG headsets with just 2 electrodes which makes it ideal to use as a wake action for AR/VR glasses and headsets. 
* The blink duration parameter in the algorithm can be used to specify different actions based on the blink duration. 