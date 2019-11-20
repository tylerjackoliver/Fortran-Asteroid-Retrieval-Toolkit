% Loop through tester.txt - get the target states of interest

function loopy()

    load reqIdxSorted.mat reqIdxSorted
    fId = fopen('tester_3435539.dat', 'r');
    fIdW = fopen('test_conditions_3435539.txt', 'w');
    i = 1;
    toCheck = 1;

    while( ~feof(fId) )

        dum = fgets(fId);

        if (i == reqIdxSorted(toCheck))

            fprintf(fIdW, '%s', dum);
            toCheck = toCheck+1;

        end

        fprintf("Completed %d\n", i);

        i = i + 1;

    end

    fclose(fId);
    fclose(fIdW);
    
end