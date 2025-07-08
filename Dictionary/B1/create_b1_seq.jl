using KomaMRI, Printf

seq = read_seq("mpf_001_new_short_copy.seq")

function scale_rf_amplitudes(input_path::String, output_path::String, scale::Float64)

    lines = readlines(input_path)
    inside_rf = false
    output_lines = String[]

    for line in lines

        if line == "[RF]"
            inside_rf = true
            push!(output_lines, line)
            println("Inside")
            continue
        end

        if inside_rf
            if isempty(line)
                inside_rf = false
                push!(output_lines, line)
                println("Outside")
                continue
            end

            parts = split(line)
            amplitude = parse(Float64, parts[2])
            parts[2] = @sprintf("%.10f", amplitude * scale)
            push!(output_lines, join(parts, " "))
            continue
        end

        push!(output_lines, line)
    end

    write(output_path, join(output_lines, "\n"))
    println("Done")
end

scale_rf_amplitudes("mpf_001_new_short_copy.seq", "mpf_001_b1scale_0.9.seq", 0.9)
