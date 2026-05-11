using JLD2
using Plots
"""
 Open QSP checkpoint file and display contents
"""
function open_qsp_chk()
    data_path = "/Users/admin/chkf_checkpoint.jld2"

    jldopen(data_path, "r") do file
        out = file["out"]

        ReCt = out["ReCt"]
        ImCt = out["ImCt"]
        intervals = out["intervals"]
        norms = out["norms"]

        println(ReCt)
        println(intervals)

        plot(intervals, ReCt,
            label="ReCt",
            xlabel="time",
            ylabel="C(t)",
            lw=2)

        #plot!(intervals, ImCt,
        #    label="ImCt",
        #    lw=2)
    end
end


open_qsp_chk()