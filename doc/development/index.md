# Developers documentation

## Installing megSAP for development (on Linux or Wiundows WSL)

In order to run the megSAP function and tool tests, perform the following steps:

- Default installation as described [here](../install_unix.md)
- execute `php src/Tools/data_setup.php -build GRCh38 -include_container` to set up data in `local_data` folder
- set up a NGSD test instance and add the credentials to the `settings.ini` file:
	
	> sudo apt-get install mariadb-server  
	> sudo mysql  
	mysql> CREATE USER 'ngsd_test'@'%' IDENTIFIED BY 'password';  
	mysql> CREATE DATABASE ngsd_test;  
	mysql> GRANT ALL ON ngsd_test.* TO 'ngsd_test'@'%';  