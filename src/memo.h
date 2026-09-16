/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_MEMO_H_
#define CARTA_ZARR_SRC_MEMO_H_

#include <map>
#include <mutex>
#include <optional>
#include <utility>

namespace carta::zarr::internal {

/**
 * A table of values computed on first use and remembered thereafter.
 *
 * The lock is held across the computation, so a value is computed exactly once no matter how many
 * threads ask for it at the same time. That also means a computation must not reach back into the
 * same table: a caller whose computation consults another Memo has to order those tables and keep
 * every caller to the same order.
 *
 * A failed computation is remembered like any other value. That is deliberate for a read-only view
 * of a store: a node that was missing stays missing for the life of the view, so what a caller
 * observes does not change under it.
 *
 * Values are handed back by reference. An entry is written once and never erased or replaced, and a
 * std::map does not move the ones it already holds when another is inserted, so a reference stays
 * valid for as long as the table does -- which is what makes a table of parsed JSON documents worth
 * having at all. A caller that stores one keeps the table alive for at least as long.
 */
template <typename Key, typename Value>
class Memo {
public:
    Memo() = default;
    Memo(const Memo&) = delete;
    Memo& operator=(const Memo&) = delete;
    Memo(Memo&&) = delete;
    Memo& operator=(Memo&&) = delete;
    ~Memo() = default;

    template <typename Compute>
    const Value& GetOrCompute(const Key& key, Compute compute) const {
        std::scoped_lock const lock(_mutex);
        const auto found = _entries.find(key);
        if (found != _entries.end()) {
            return found->second;
        }
        const auto insertion = _entries.emplace(key, compute());
        return insertion.first->second;
    }

    // Put a value in that was not computed here. It is for a table whose entries are already in
    // hand -- the node metadata a store's consolidated copy arrives holding -- so that they are one
    // value looked up one way rather than a second copy consulted first. An entry that is already
    // there wins, because it is the one callers may be holding a reference to.
    void Insert(Key key, Value value) const {
        std::scoped_lock const lock(_mutex);
        _entries.emplace(std::move(key), std::move(value));
    }

private:
    mutable std::mutex _mutex;
    mutable std::map<Key, Value> _entries;
};

/**
 * One value computed on first use and remembered thereafter: Memo with nothing to key on.
 *
 * Shares Memo's contract, including that a failed computation is remembered.
 */
template <typename Value>
class Lazy {
public:
    Lazy() = default;
    Lazy(const Lazy&) = delete;
    Lazy& operator=(const Lazy&) = delete;
    Lazy(Lazy&&) = delete;
    Lazy& operator=(Lazy&&) = delete;
    ~Lazy() = default;

    template <typename Compute>
    const Value& GetOrCompute(Compute compute) const {
        std::scoped_lock const lock(_mutex);
        if (!_value.has_value()) {
            _value.emplace(compute());
        }
        return *_value;
    }

private:
    mutable std::mutex _mutex;
    mutable std::optional<Value> _value;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_MEMO_H_
