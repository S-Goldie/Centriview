## Instructions for Streamlit Deployment

_These instructions are only for users wishing to deploy their own local version of the Centriview WebApp_

* Download all .py and .png files, maintaining folder structure. 
* Deploy the Home.py file using `streamlit run Home.py` in a python console.
* A browser window should open automatically. The local and network address will be shown in the Python console.
* If the correct firewall port is opened on the local machine (8501 by default), other users on the local network should be able to use the local network IP address to access the webapp through standard browsers.
* The browser can be closed, but if the local Python environment running streamlit is closed then naturally the webapp will close for all users.
