
CPLEX_studio_binaries_key="CPLEX_STUDIO_BINARIES"

if haskey(ENV, CPLEX_studio_binaries_key)
    if isdir(ENV[CPLEX_studio_binaries_key])
        @info "Checked 1 location, found 1 location for CPLEX, using $(ENV[CPLEX_studio_binaries_key])."
        println(ENV[CPLEX_studio_binaries_key])
        exit()
    end
end
locations = Dict(
    "~/.local/bin/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/bin/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/bin/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/bin/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/bin/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/bin/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "~/.local/share/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/share/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/share/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/share/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/share/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/share/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "~/.local/lib/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/lib/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/lib/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/lib/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "~/.local/lib/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "~/.local/lib/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "/usr/lib/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/lib/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/lib/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/lib/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/lib/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/lib/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "/usr/share/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/share/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/share/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/share/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/share/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/share/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "/usr/local/lib/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/lib/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/local/lib/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/lib/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/local/lib/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/lib/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "/usr/local/share/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/share/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/local/share/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/share/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/usr/local/share/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/usr/local/share/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    # ---
    "/opt/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/opt/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/opt/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/opt/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
    "/opt/ibm/CPLEX_Studio201/cplex/bin/x86-64_linux/" => false,
    "/opt/ibm/CPLEX_Studio1210/cplex/bin/x86-64_linux/" => false,
)
@assert !isempty(locations) "Locations should not be empty before the search."
prel = length(locations)


for (location, cplexpresent) in locations
    locations[location] = cplexpresent || isdir(location)
end

filter!(l -> l.second, locations)

if isempty(locations)
    # try to find other locations
    error("No CPLEX installation found.")
end
@assert !isempty(locations) "Locations should not be empty."

l = length(locations)
loc = first(locations).first

if l > 1
    # Do some special handling for multiple files.
end

@info "Checked $prel location$(prel == 1 ? "" : "s"), found $l location$(l == 1 ? "" : "s") for CPLEX, using \"$loc\"."

println(loc)
