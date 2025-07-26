using Printf

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

seq = "C:/Users/lucya/MSC_PROJECT/Dictionary/sequences/mpf_001_PhantomStudy_short_124.seq"
b1s = vcat(0.8:0.02:0.9, 0.91:0.01:1.09, 1.10:0.02:1.2)
for b1 in b1s
    filename = @sprintf("b1_%.2f_short_124.seq", b1)
    scale_rf_amplitudes(seq, filename, b1)
end
