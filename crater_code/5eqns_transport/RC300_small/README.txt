
Test problem for smaller domain and short times.

To run:

    make # to create xclaw

    mkdir _output   # or desired name for this run
    cp *.data _output
    cd _output

    # modify _output/claw2ez.data to modify resolution by changing mx,my
    # change method(9) to switch from muscl (1) to thinc (2).

    time ../xclaw  # to run and report timing

To plot:

    Need to clean up Python plotting script
