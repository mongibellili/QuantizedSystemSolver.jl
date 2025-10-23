






function push_scope!(s::SymbolTableStack, table::SymbolTable)
    push!(s.tables, table)
end
function pop_scope!(s::SymbolTableStack)
    isempty(s.tables) && error("No more scopes to pop")
    pop!(s.tables)
end
function peek(s::SymbolTableStack)
    isempty(s.tables) && error("No active scope to peek into")
    return s.tables[end]
end


# Lookup a symbol starting from currentTable (innermost scope) to bottom (outermost)
function lookup(s::SymbolTableStack, name::Symbol)
    for table in Iterators.reverse(s.tables)
        if haskey(table.entries, name)
            return table.entries[name]
        end
    end
    return nothing  # Not found in model or user defined as global
end

#safe_to_inline is need in the case of a local unsafe_to_swap param and a global safe_to_inline; we want to store this unsafe param in table so that a local equation sees this and does not swap the global
function add_symbol!(table::SymbolTable, name::Symbol, value, safe_to_inline::Bool=true)
    value= safe_to_inline ? value : :nothing # If not safe to inline, use :nothing. This prevents inlining but still tracks the symbol.
    entry = SymbolEntry(value, safe_to_inline)
    table.entries[name] = entry
end










