<!-- #region -->
## Running the Course Image with Singularity
Docker requires root permissions to run, so you are unlikely to be able to run Docker on a computer that you are not fully in control of.  As an alternative you can run the course image with [Singularity](https://sylabs.io/singularity/), another container system. Singularity is similar to Docker, and can run Docker images, but you do not need special permissions to run Singularity images *or* Docker images with Singularity (as long as Singularity is actually installed on the computer).


### Running the Course Image on a SLURM cluster

We will use the example of the Duke Computer Cluster, but these instructions should be easily adaptable to other clusters

1. From your computer run this to connect to DCC:
```
ssh NetID@dcc-login-03.oit.duke.edu
```
2. Once you are connected run this to start a tmux session:
```
tmux new -s jupyter
```
3. Once you have started a tmux session you can start up Jupyter by running this command in the *top directory of this repository*:
```
srun -A chsi -p chsi --cpus-per-task=5 --mem=10G admin/run_singularity/run_singularity_jupyter.sh
```
> Note: the first time you run this, it might take a VERY long time to download the Docker image and build the Singularity image from it

Running this command will print a bunch of stuff. You can ignore everything except the last two lines, which will say something like:
```
http://dcc-chsi-01:8889/?token=08172007896ad29bb5fbd92f6f3f516a8b2f7303ed7f1df3
or http://127.0.0.1:8889/?token=08172007896ad29bb5fbd92f6f3f516a8b2f7303ed7f1df3
```

You need this information for the next few steps. For the next step you need the “dcc-chsi-01:8889” part.
“dcc-chsi-01” is the compute node that Jupyter is running on and “8889” is the port it is listening on. You may get a different value every time you start the container.

4. You want to run the following command in another terminal on your computer to set up port forwarding.
```
ssh -L PORT:NODE.rc.duke.edu:PORT NetID@dcc-login-03.oit.duke.edu
```
In this command you want to replace “PORT” with the value you got for port from the srun command and replace “NODE” with the compute node that was printed by the srun command. So for the example above, the ssh port forwarding command would be:

```
ssh -L 8889:dcc-chsi-01.rc.duke.edu:8889 NetID@dcc-login-03.oit.duke.edu
```

5. Now you can put the last line that the srun command printed in your web browser and it should open your Jupyter instance running on DCC.

#### Notes
1. The Jupyter session keeps running until you explicitly shut it down.  If the port forwarding SSH connection drops you will need to restart SSH with the same command, but you don’t need to restart Jupyter.

2. There are two ways to explicitly shut down Jupyter:
    1. Within Jupyter, click on the *Jupyter* logo in the top left to go to the main Jupyter page, then click "Quit" in the top right
    2. Do control-C twice in the terminal where you started Jupyter. If this connection dropped, you can reconnect to it with:
    ```
    ssh NetID@dcc-login-03.oit.duke.edu
    tmux a -t jupyter
    ```
    After shutting down the Jupyter session you can type `exit` at the terminal to close the tmux session.

3. If you need more memory or more cpus you can use the `--mem` and/or `--cpus-per-task` arguments to in the “srun”, for example to request 4 CPUs and 10GB of RAM:
```
srun --cpus-per-task=4 --mem=10G admin/run_singularity/run_singularity_jupyter.sh
```

4. If you have high priority access to a partition you can request that partition be used with the `-A` and `-p` arguments to `srun`:
```
srun -A chsi -p chsi admin/run_singularity/run_singularity_jupyter.sh
```

6. You can combine several of these command line flags:
```
srun -A chsi -p chsi --cpus-per-task=4 --mem=10G admin/run_singularity/run_singularity_jupyter.sh
```

6. It is strongly recommended to set the `SINGULARITY_CACHEDIR` environment variables in your .bashrc or when running `srun`. This environment variable specifies where the Docker image (and the Singularity image built from it) are saved. If this variable is not specified, singularity will cache images in `$HOME/.singularity/cache`, which can fill up quickly. This is discussed in the [Singularity Documentation](https://sylabs.io/guides/3.7/user-guide/build_env.html#cache-folders)

```
mkdir -p "/work/${USER}"
export SINGULARITY_CACHEDIR="/work/${USER}/singularity_cache"; srun -A chsi -p chsi --cpus-per-task=4 --mem=10G admin/run_singularity/run_singularity_jupyter.sh
```

### Install Singularity
Here are instructions for installing:

- [Singularity version 2.6](https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps)
- [Singularity version 3.2](https://sylabs.io/guides/3.2/user-guide/quick_start.html#quick-installation-steps)
- [Singularity Desktop for macOS (Alpha Preview)](https://sylabs.io/singularity-desktop-macos/)
<!-- #endregion -->
